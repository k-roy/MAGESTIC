
library(tidyverse)
designed_oligos <- read_tsv('/path/to/scripts_and_keyfiles/by_project/NNS/bc1_donor_bc0/keyfiles/20240411_Twist_200mer_oligo_array_order.tsv')

df <- read_tsv('20210906_and_20210820_SpCas9_LbCas12a_QTL_plasmid_libraries_bc0_counts_guide_donor_matches_minimal_columns.tsv')

df %>% group_by(scaffold) %>% tally() %>% ggplot(aes(x = V_library))

df <- df %>% 
mutate(sequence_perfect = perfect_leader & perfect_guide & perfect_donor,
  sequence_perfect_mid_donor = perfect_leader & perfect_guide & perfect_middle_donor_region) 

ggplot(df, aes(x = V_library)) + geom_bar() + facet_wrap(sequence_perfect ~ sequence_perfect_mid_donor)

df_3_cts_min <- df %>% filter(counts > 3)

df_3_cts_min_summary <- df_3_cts_min %>% group_by(scaffold) %>% tally()

perfect_guide_donor_summary <- df_3_cts_min %>% 
  filter(perfect_leader, perfect_guide, perfect_donor) %>% 
  group_by(scaffold) %>% tally()

perfect_middle_donor_region_summary <- df_3_cts_min %>% 
  filter(perfect_leader, perfect_guide, perfect_middle_donor_region) %>% 
  group_by(scaffold) %>% tally()

df_3_cts_min_summary <- df_3_cts_min_summary %>% 
  left_join( perfect_middle_donor_region_summary %>% rename(perfect_mid_donor = n) ) %>%
  left_join( perfect_guide_donor_summary %>% rename(perfect_guide_donor = n) ) %>%
  rename(total = n)
  
df_3_cts_min_summary <- df_3_cts_min_summary %>% mutate(fraction_perfect = perfect_guide_donor / total, 
                                fraction_middle_donor_perfect = perfect_mid_donor / total)

dput(as.character(unique(df_3_cts_min_summary$scaffold)))

df_3_cts_min_summary$scaffold <- factor(df_3_cts_min_summary$scaffold,
                                        levels = c("WT_SpCas9",
                                                   "WT_LbCas12a",
                                                   "GC1_LbCas12a", 
                                                   "GC13_LbCas12a", 
                                                   "mut_LbCas12a" 
                                                   ))

pdf('/labs/mylab/kevinroy/processed_data/20210822_SpCas9_LbCas12a_QTL_plasmid_libraries/plots/step_1_libraries_3_counts_min_summary_by_scaffold.pdf', width = 6, height = 6)
ggplot(df_3_cts_min_summary, aes(x = scaffold, y = total)) + geom_col() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Number of unique barcodes across step 1 libraries\n(min. 3 reads for bc0 guide-donors) ')
dev.off()

df_3_cts_min_summary_long <- df_3_cts_min_summary %>% 
  pivot_longer(cols = c("fraction_perfect", "fraction_middle_donor_perfect"), 
               names_to = "category") %>% rename(fraction = value)

df_3_cts_min_summary_long$category = factor(df_3_cts_min_summary_long$category,
                                            levels = c("fraction_perfect", "fraction_middle_donor_perfect"))

pdf('/labs/mylab/kevinroy/processed_data/20210822_SpCas9_LbCas12a_QTL_plasmid_libraries/plots/step_1_libraries_3_cts_min_summary_fraction_perfect.pdf', width = 6, height = 6)

ggplot(df_3_cts_min_summary_long %>% filter(scaffold != "mut_LbCas12a" ), 
       aes(x = scaffold, y = fraction, fill = category)) + 
  geom_col(position = position_dodge()) + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Fraction sequence perfect guide-donors\n(min. 3 reads for unique bc0 guide-donors) ')
dev.off()

oligo_info <- read_tsv('SpCas9_LbCas12a_impLbCas12a_QTL_oligos_info.txt', guess_max = 10000000)

oligo_info %>% group_by(codon_variant_type) %>% tally() %>% print(n = 1000)

df_3_cts_min_info <- df %>% left_join(oligo_info, by = c("oligo_name"))

ann_oligo_info <- oligo_info %>% mutate(in_library = oligo_name %in% df_3_cts_min_info$oligo_name)
ann_oligo_info$in_library

ann_oligo_info %>% filter(!in_library) %>% 
  write_tsv('guide_donors_not_detected_in_step_1_libraries.tsv')

df_3_cts_min_info %>% filter(is.na(codon_variant_type)) 


# df_3_cts_min_info %>% group_by(V_library, scaffold) %>% select(oligo_name) %>% unique() %>% tally()

df_3_cts_min_info %>% group_by(scaffold, codon_variant_type) %>% tally() %>% print(n = 1000)

scaffold_by_variant_type <- df_3_cts_min_info %>% filter(perfect_leader, perfect_guide, perfect_donor) %>%
  group_by(scaffold) %>% select(oligo_name) %>% unique() %>% tally()

unique(ann_oligo_info$category)
dput(as.character(unique(ann_oligo_info$category)))

pdf('/labs/mylab/kevinroy/processed_data/20210822_SpCas9_LbCas12a_QTL_plasmid_libraries/plots/fraction_library_covered_3_cts_min_summary.pdf', width = 6, height = 6)
ggplot(ann_oligo_info %>% filter(category %in% c("codon_changes", "natural_variants")), aes( x = codon_variant_type, fill = in_library)) +
  geom_bar() + facet_grid(nuclease_PAM ~ category) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Fraction library covered\n(min. 3 reads for unique bc0 guide-donors) ')
dev.off()

ggplot(ann_oligo_info %>% filter(category %in% c("PDR5_controls", "DHY214_controls")), 
       aes( x = codon_variant_type, fill = in_library)) +
  geom_bar() + facet_grid(nuclease_PAM ~ category) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Fraction library covered\n(min. 3 reads for unique bc0 guide-donors) ')

ggplot(ann_oligo_info %>% filter(category %in% c("PTC_del_syn_controls")), 
       aes( x = codon_variant_type, fill = in_library)) +
  geom_bar() + facet_grid(nuclease_PAM ~ category) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Fraction library covered\n(min. 1 read for unique bc0 guide-donors) ')
