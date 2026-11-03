# NRD1 fitness screen: fitness index and hit calls (Supplementary Table 1, Figure 1D)

```bash
magestic-pooled-slope --sets sets.json --reference <bc1 reference table> --config config.json --out-dir out/
```

`sets.json` is `sets.template.json` with `<SETS>` pointing at the six per-set normalized count tables
(CPM per set, `Barcode, Count_<g>, Normalized_Count_<g>`). The consolidation adds a pseudocount of 1 to every
normalized count, fills a barcode absent from a whole set with 0, and keeps only barcodes that resolve to a single
donor. `config.json` holds the settings Supplementary Table 1 used: `tail_prob` 0.05 and `depletion_fc` 0.1, where
the original R script (`NRD1/fitness_slope_heatmap.R`) defaults to 0.02 and 1e-4.

`hits_with_support.csv` holds the edits that meet criteria 1 to 3 and have at least one barcode meeting them on
its own (criterion 4). The other files reproduce the R script's tables under its file names.

Measured (planning/380): on the consolidated input the original analysis used, and on the two inputs of
planning/379, all six tables equal the R script's output to within 1.8e-15, and the hit set is the 23 edits of
Supplementary Table 1.
