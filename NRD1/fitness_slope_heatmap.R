## ============================== SETUP ===================================== ##
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(rlang)
})

# -------- CONFIG --------
setwd("~/Desktop/Bc1 analysis/")

in_file  <- "consolidated_normalized_counts_YPD.csv"  # input CSV
out_dir  <- "."                                        # output folder

# analysis knobs
min_total_count  <- 10     # keep barcode if avg G0 >= this
eps              <- 1e-8   # pseudocount to avoid log(0)
use_weighted     <- TRUE   # weight by G0 when collapsing across barcodes
fc_trim_cutoff   <- 0.01   # drop series after it hits k consecutive points below this FC
sustained_k      <- 3      # number of consecutive points below cutoff (post-G0) to trigger a trim

# hit calling knobs 
treat_ypd        <- "YPD"
treat_ypd_noIA   <- "YPD_noIAA"
hit_tail_prob    <- 0.02   # 2nd percentile cutoff for YPD adjusted slopes
diff_factor_rule <- 3      # require slope_adj_YPD / 3 < slope_adj_YPD_noIAA
depletion_fc_cut <- 1e-4   # depletion in last FC for YPD

# plot outputs
heatmap_slopes_pdf <- file.path(out_dir, "heatmap_log2FC_slopes.pdf")
heatmap_time_pdf   <- file.path(out_dir, "heatmap_log2FC_over_time.pdf")
dispersion_pdf     <- file.path(out_dir, "dispersion_log2FC_slopes.pdf")
regression_pdf     <- file.path(out_dir, "regression_line.pdf")

# tables
slope_table_csv    <- file.path(out_dir, "slope_table_log2FC.csv")
hits_csv           <- file.path(out_dir, "hits_summary.csv")
mean_time_csv      <- file.path(out_dir, "mean_log2FC_by_group.csv")
depletion_csv      <- file.path(out_dir, "depletion_log2FC.csv")
bc_flags_csv       <- file.path(out_dir, "barcode_support_flags.csv")
support_counts_csv <- file.path(out_dir, "hit_barcode_support_summary.csv")

## ============================ LOAD & CHECKS =============================== ##
wide <- readr::read_csv(in_file, show_col_types = FALSE) %>%
  mutate(
    bc1 = as.character(bc1),
    amino_acid_position = suppressWarnings(as.integer(amino_acid_position))
  )

# keep only bc0 fragments passing filter
if ("multiple_donor_bc0_fragments_passing_filter" %in% names(wide)) {
  wide <- wide %>% filter(.data$multiple_donor_bc0_fragments_passing_filter == "TRUE")
}

# ---- mutation classification (PTC / synonymous / missense / NA) ----
trim_chr <- function(x) if (is.character(x)) trimws(x) else x
ref <- trim_chr(wide$ref_amino_acid)
mut <- trim_chr(wide$mutant_amino_acid)
mutation_type <- ifelse(is.na(ref) | is.na(mut), NA_character_,
                        ifelse(mut == "*", "PTC",
                               ifelse(ref == mut, "synonymous", "missense")))

# insert mutation_type right before the first *_Normalized_Count_* column, if present
norm_cols  <- grep("_Normalized_Count_", names(wide), fixed = TRUE)
insert_pos <- if (length(norm_cols) > 0) min(norm_cols) else ncol(wide) + 1
wide <- if (insert_pos == 1) {
  cbind(mutation_type = mutation_type, wide)
} else if (insert_pos > ncol(wide)) {
  cbind(wide, mutation_type = mutation_type)
} else {
  cbind(wide[, 1:(insert_pos - 1), drop = FALSE],
        mutation_type = mutation_type,
        wide[, insert_pos:ncol(wide), drop = FALSE])
}

## ======================== LONG FORMAT & PARSING =========================== ##
long <- wide %>%
  pivot_longer(cols = tidyselect::matches("_Normalized_Count_"),
               names_to = "sample", values_to = "fitness") %>%
  mutate(
    treatment_raw = sub("_Normalized_Count_.*", "", sample),
    replicate     = str_extract(treatment_raw, "_R\\d+$") %>% str_remove("_R") %>% as.integer(),
    set_flag      = str_detect(treatment_raw, "_set$"),
    treatment     = treatment_raw %>%
      sub("_R\\d+$", "", .) %>%
      sub("_set$", "", .) %>%
      gsub("\\s+", "_", .),
    generation    = as.integer(str_extract(sample, "(?<=Normalized_Count_)\\d+"))
  ) %>%
  select(bc1, gene_name, amino_acid_position, mutation_type,
         sample, treatment, set_flag, replicate, generation, fitness) %>%
  filter(!is.na(amino_acid_position))

## ===================== PER-BARCODE NORMALIZATION ========================== ##
# average G0 per barcode (across replicates within bc1)
gen0_means <- long %>%
  filter(generation == 0) %>%
  group_by(bc1) %>%
  summarise(g0 = mean(fitness, na.rm = TRUE), .groups = "drop")

valid_bc <- gen0_means %>%
  filter(is.finite(g0), g0 >= min_total_count) %>%
  pull(bc1)

df <- long %>%
  filter(bc1 %in% valid_bc) %>%
  left_join(gen0_means, by = "bc1") %>%
  mutate(
    # harmonize gen0 to the per-barcode average
    fitness = if_else(generation == 0, g0, fitness),
    fc      = pmax(fitness, eps) / pmax(g0, eps),
    log2fc  = log2(fc),
    w       = if (use_weighted) pmax(g0, 1) else 1
  )

## ============================ HELPERS ===================================== ##
# Single trimming rule: find first run of k consecutive post-G0 points with log2FC < log2(cutoff_fc);
# keep up to and including that generation. NA does NOT count toward the run.
trim_series <- function(g, y, cutoff_fc = 0.5, k = 2) {
  cutoff_l2 <- log2(cutoff_fc)
  ord <- order(g); g <- g[ord]; y <- y[ord]
  finite <- is.finite(y) & is.finite(g)
  g <- g[finite]; y <- y[finite]
  if (length(g) == 0) return(tibble(generation = numeric(0), y = numeric(0)))
  if (!any(g == 0, na.rm = TRUE)) return(tibble(generation = g, y = y))  # no G0 → no trim
  
  post_mask <- g > 0
  gens_post <- g[post_mask]
  y_post    <- y[post_mask]
  if (length(y_post) == 0) return(tibble(generation = g, y = y))
  
  below <- y_post < cutoff_l2
  below_flag <- ifelse(is.na(below), FALSE, below)
  cnt <- 0L; first_idx <- NA_integer_
  for (i in seq_along(below_flag)) {
    cnt <- if (below_flag[i]) cnt + 1L else 0L
    if (cnt >= k) { first_idx <- i; break }
  }
  if (is.na(first_idx)) return(tibble(generation = g, y = y))
  cut_gen <- gens_post[first_idx]
  keep <- g <= cut_gen
  tibble(generation = g[keep], y = y[keep])
}

safe_lm_slope <- function(x, y, min_points = 3) {
  d <- tibble(x = as.numeric(x), y = as.numeric(y)) %>%
    filter(is.finite(x), is.finite(y)) %>%
    distinct(x, .keep_all = TRUE) %>% arrange(x)
  if (nrow(d) < min_points) return(NA_real_)
  coef(lm(y ~ x, data = d))[2]
}

## ===================== COLLAPSE (MEAN) & TRIM ONCE ======================== ##
# Collapsed mean log2FC by variant×treatment×generation
all_summaries_fc <- df %>%
  group_by(amino_acid_position, mutation_type, treatment, generation) %>%
  summarise(
    mean_log2fc = if (use_weighted) weighted.mean(log2fc, w, na.rm = TRUE)
    else mean(log2fc, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(all_summaries_fc, mean_time_csv)

# Determine cut_gen from the COLLAPSED series (one cut per variant×treatment)
trim_info <- all_summaries_fc %>%
  group_by(amino_acid_position, mutation_type, treatment) %>%
  summarise(
    cut_gen = {
      ser <- trim_series(generation, mean_log2fc,
                         cutoff_fc = fc_trim_cutoff, k = sustained_k)
      if (nrow(ser) == 0) max(generation[is.finite(mean_log2fc)], na.rm = TRUE)
      else max(ser$generation, na.rm = TRUE)
    },
    .groups = "drop"
  )

# Fit slopes ONLY on collapsed, trimmed series
collapsed_slopes <- all_summaries_fc %>%
  inner_join(trim_info, by = c("amino_acid_position","mutation_type","treatment")) %>%
  filter(is.finite(mean_log2fc), is.finite(generation), generation <= cut_gen) %>%
  group_by(amino_acid_position, mutation_type, treatment) %>%
  summarise(
    slope = safe_lm_slope(generation, mean_log2fc),
    .groups = "drop"
  )

# Median-center within treatment
slope_table <- collapsed_slopes %>%
  group_by(treatment) %>%
  mutate(
    median_slope = median(slope, na.rm = TRUE),
    slope_adj    = slope - median_slope
  ) %>%
  ungroup()

## ============================ DEPLETION =================================== ##
# For YPD only: compare gen0 vs last finite mean point
last_by_group <- all_summaries_fc %>%
  group_by(amino_acid_position, mutation_type, treatment) %>%
  summarise(last_gen = max(generation[is.finite(mean_log2fc)], na.rm = TRUE),
            .groups = "drop")

fitness_extremes <- all_summaries_fc %>%
  filter(generation == 0) %>%
  select(amino_acid_position, mutation_type, treatment, log2fc_gen0 = mean_log2fc) %>%
  right_join(
    all_summaries_fc %>%
      inner_join(last_by_group, by = c("amino_acid_position","mutation_type","treatment")) %>%
      filter(generation == last_gen) %>%
      select(amino_acid_position, mutation_type, treatment, log2fc_last = mean_log2fc),
    by = c("amino_acid_position","mutation_type","treatment")
  ) %>%
  mutate(depletion = is.finite(log2fc_last) & (log2fc_last < log2(depletion_fc_cut))) %>%
  select(amino_acid_position, mutation_type, treatment, log2fc_gen0, log2fc_last, depletion)

readr::write_csv(fitness_extremes, depletion_csv)

## ============================= HIT CALLING ================================= ##
# Wide compare for YPD vs YPD_noIAA
slope_wide <- slope_table %>%
  filter(treatment %in% c(treat_ypd, treat_ypd_noIA),
         mutation_type %in% c("missense", "synonymous", "PTC", "Deletion")) %>%
  select(amino_acid_position, mutation_type, treatment, slope, median_slope, slope_adj) %>%
  pivot_wider(names_from = treatment, values_from = c(slope, median_slope, slope_adj), names_sep = "_")

# bring depletion (YPD only)
depl_by_posmut <- fitness_extremes %>%
  filter(treatment == treat_ypd) %>%
  select(amino_acid_position, mutation_type, depletion)

slope_wide <- slope_wide %>%
  left_join(depl_by_posmut, by = c("amino_acid_position","mutation_type"))

# threshold from YPD distribution (adjusted slopes)
hit_threshold_ypd <- quantile(
  slope_table$slope_adj[slope_table$treatment == treat_ypd],
  probs = hit_tail_prob, na.rm = TRUE
)

hits <- slope_wide %>%
  filter(!is.na(slope_adj_YPD), !is.na(slope_adj_YPD_noIAA), depletion == TRUE) %>%
  filter(
    slope_adj_YPD < hit_threshold_ypd,
    slope_adj_YPD / diff_factor_rule < slope_adj_YPD_noIAA
  ) %>%
  mutate(hit = TRUE)

## ================= BARCODE-LEVEL SUPPORT ========= ##
# Use the SAME cut_gen window (from collapsed series) for all barcodes in the group
per_bc_slopes <- df %>%
  inner_join(trim_info, by = c("amino_acid_position","mutation_type","treatment")) %>%
  filter(is.finite(log2fc), is.finite(generation), generation <= cut_gen) %>%
  group_by(amino_acid_position, mutation_type, treatment, bc1) %>%
  summarise(
    slope = safe_lm_slope(generation, log2fc),
    w = first(w),
    .groups = "drop"
  )

# Per-barcode adjusted slopes (median-center by treatment)
median_by_treat_bc <- per_bc_slopes %>%
  group_by(treatment) %>%
  summarise(median_slope_bc = median(slope, na.rm = TRUE), .groups = "drop")

per_bc_slopes_adj <- per_bc_slopes %>%
  left_join(median_by_treat_bc, by = "treatment") %>%
  mutate(slope_adj = slope - median_slope_bc)

# Per-barcode depletion (YPD only), using the last finite log2FC for that barcode
bc_depl_ypd <- df %>%
  filter(treatment == treat_ypd) %>%
  arrange(amino_acid_position, mutation_type, bc1, generation) %>%
  group_by(amino_acid_position, mutation_type, bc1) %>%
  summarise(
    log2fc_last = {
      lf <- log2fc[is.finite(log2fc)]
      gg <- generation[is.finite(log2fc)]
      if (length(lf) == 0) NA_real_ else lf[which.max(gg)]
    },
    depletion_bc = is.finite(log2fc_last) & (log2fc_last < log2(depletion_fc_cut)),
    .groups = "drop"
  )

# Wide per-barcode adjusted slopes for YPD vs YPD_noIAA (no barcode-specific trimming)
bc_slope_wide <- per_bc_slopes_adj %>%
  filter(treatment %in% c(treat_ypd, treat_ypd_noIA)) %>%
  select(amino_acid_position, mutation_type, bc1, treatment, slope_adj) %>%
  pivot_wider(names_from = treatment, values_from = slope_adj) %>%
  left_join(bc_depl_ypd, by = c("amino_acid_position", "mutation_type", "bc1"))

# dynamic column name handling (in case pivot names differ)
col_ypd   <- if ("slope_adj_YPD" %in% names(bc_slope_wide)) "slope_adj_YPD" else treat_ypd
col_noIAA <- if ("slope_adj_YPD_noIAA" %in% names(bc_slope_wide)) "slope_adj_YPD_noIAA" else treat_ypd_noIA
stopifnot(col_ypd %in% names(bc_slope_wide), col_noIAA %in% names(bc_slope_wide))

bc_hit_flags <- bc_slope_wide %>%
  mutate(
    supports_hit = !is.na(.data[[col_ypd]]) &
      !is.na(.data[[col_noIAA]]) &
      depletion_bc &
      (.data[[col_ypd]] < hit_threshold_ypd) &
      (.data[[col_ypd]] / diff_factor_rule < .data[[col_noIAA]]),
    agrees_direction = !is.na(.data[[col_ypd]]) &
      !is.na(.data[[col_noIAA]]) &
      (.data[[col_ypd]] < .data[[col_noIAA]]) &
      (.data[[col_ypd]] < 0)
  )

support_counts <- bc_hit_flags %>%
  group_by(amino_acid_position, mutation_type) %>%
  summarise(
    n_barcodes  = dplyr::n(),
    n_support   = sum(supports_hit, na.rm = TRUE),
    n_agree     = sum(agrees_direction, na.rm = TRUE),
    pct_support = ifelse(n_barcodes > 0, n_support / n_barcodes, NA_real_),
    pct_agree   = ifelse(n_barcodes > 0, n_agree   / n_barcodes, NA_real_),
    .groups = "drop"
  )

readr::write_csv(bc_hit_flags, bc_flags_csv)
readr::write_csv(support_counts, support_counts_csv)

## ================= MERGE & SAVE TABLES ==================================== ##
slope_table <- slope_table %>%
  left_join(support_counts, by = c("amino_acid_position", "mutation_type"))

hits <- hits %>%
  left_join(support_counts, by = c("amino_acid_position", "mutation_type"))

# Add hit flag back to slope_table
slope_table <- slope_table %>%
  left_join(hits %>% select(amino_acid_position, mutation_type, hit),
            by = c("amino_acid_position", "mutation_type")) %>%
  mutate(hit = if_else(is.na(hit), FALSE, hit))

# Optional pseudo "HIT" row for heatmaps
hit_rows <- slope_table %>%
  filter(hit == TRUE) %>%
  distinct(amino_acid_position, mutation_type) %>%
  mutate(treatment = "HIT", slope = NA_real_, median_slope = NA_real_, slope_adj = NA_real_, hit = TRUE)

slope_table_out <- bind_rows(slope_table, hit_rows) %>%
  mutate(treatment = factor(treatment, levels = c("HIT", treat_ypd, treat_ypd_noIA)))

# Save
readr::write_csv(slope_table_out, slope_table_csv)
readr::write_csv(hits,        hits_csv)

## ================================ PLOTS =================================== ##
# ---- Heatmap of adjusted slopes (excluding PTC/Deletion) ----
heatmap_table <- slope_table_out %>%
  filter(!mutation_type %in% c("PTC", "Deletion"))

ggplot(heatmap_table, aes(x = amino_acid_position, y = treatment, fill = slope_adj)) +
  geom_tile(color = "grey85") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "blue",
                       midpoint = 0, name = "Normalized Slope", na.value = "grey90") +
  facet_wrap(~mutation_type, ncol = 1, strip.position = "top") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Amino Acid Position", y = "Treatment",
       title = "Fitness Slope Heatmap (Median Normalized)")
ggsave(heatmap_slopes_pdf, width = 12, height = 8)

# ---- Dispersion plot: YPD_noIAA vs YPD (adjusted slopes) ----
slope_wide_plot <- slope_table %>%
  filter(treatment %in% c(treat_ypd, treat_ypd_noIA)) %>%
  select(amino_acid_position, mutation_type, treatment, slope_adj, hit) %>%
  pivot_wider(names_from = treatment, values_from = slope_adj) %>%
  drop_na(!!sym(treat_ypd), !!sym(treat_ypd_noIA))

ggplot(slope_wide_plot, aes_string(x = treat_ypd_noIA, y = treat_ypd)) +
  geom_point(aes(color = mutation_type, size = hit),
             shape = 16, alpha = 0.85, stroke = 0.5) +
  scale_color_manual(values = c("missense" = "red2", "synonymous" = "blue3",
                                "PTC" = "grey40", "Deletion" = "black")) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = paste0("Adjusted slope (", treat_ypd_noIA, ")"),
       y = paste0("Adjusted slope (", treat_ypd, ")"),
       title = "Adjusted log2FC Slopes by Mutation Type",
       color = "Mutation Type") +
  theme_minimal(base_size = 13)
ggsave(dispersion_pdf, width = 12, height = 8)