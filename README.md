# RefinedIBD LC-vs-HC evaluation

Scripts to evaluate and visualize RefinedIBD segments detected in low-coverage (`LC`) and high-coverage (`HC`) datasets.

`HC` is treated as the reference set and `LC` as the query/predicted set. The main metric script compares IBD intervals by chromosome and individual pair after removing coverage suffixes from sample IDs.

## Scripts

- `Metricas.py`: computes LC-vs-HC interval overlap metrics and plots Recall, Accuracy, and F1 by IBD segment length.
- `MosaicoIBD.R`: plots all IBD segments for one target individual, separated by IBD partner and dataset (`hcWGS` / `lcWGS`). Segments are colored by length class.
- `CircusIBD.R`: creates a circular pairwise summary plot for a target individual and one selected IBD partner.


## Running RefinedIBD

The repository includes two SLURM scripts used to run RefinedIBD by chromosome:

- `RunRefinedLC.sh`: low-coverage dataset.
- `RunRefinedHC.sh`: high-coverage dataset.

Submit them as chromosome arrays:

```bash
sbatch RunRefinedLC.sh
sbatch RunRefinedHC.sh
```

For local testing of a single chromosome, the scripts also accept a chromosome number as the first argument:

```bash
bash RunRefinedLC.sh 1
bash RunRefinedHC.sh 1
```

Each array task runs one chromosome. The scripts expect chromosome-specific VCFs and genetic maps and write outputs to `LC/` or `HC/`:

```text
LC/chr<CHR>_ibd.ibd.gz
HC/chr<CHR>_ibd.ibd.gz
```

The scripts then run `merge-ibd-segments` to reduce fragmentation of nearby RefinedIBD calls:

```text
LC/chr<CHR>_ibd_merged.ibd
HC/chr<CHR>_ibd_merged.ibd
```

The merge step uses:

```bash
java -jar /opt/merge-ibd-segments.17Jan20.102.jar VCF MAP 0.01 1 < input.ibd > output_merged.ibd
```

`merge-ibd-segments` reads the RefinedIBD `.ibd` content from standard input. Passing the input file as an extra command-line argument is incorrect and produces a `usage:` message instead of merged segments.

The merge step is included because low-coverage imputed data can fragment an underlying IBD tract into multiple adjacent calls. Merging helps reduce this fragmentation, but it does not solve the main biological/technical limitation observed here: LC still produces many false positives. A plausible explanation is that imputation can reconstruct the same missing or low-information regions similarly across unrelated individuals, creating apparent shared haplotypes by chance. In that situation, merging adjacent fragments may make the calls cleaner but cannot distinguish true IBD from false sharing introduced by imputation artifacts.

## Input files

The RefinedIBD SLURM scripts first generate the standard compressed outputs:

```text
HC/chr1_ibd.ibd.gz
HC/chr2_ibd.ibd.gz
...
LC/chr1_ibd.ibd.gz
LC/chr2_ibd.ibd.gz
...
```

They then generate merged outputs used by the final metric and visualization scripts:

```text
HC/chr1_ibd_merged.ibd
HC/chr2_ibd_merged.ibd
...
LC/chr1_ibd_merged.ibd
LC/chr2_ibd_merged.ibd
...
```

Each IBD file is expected to follow the RefinedIBD-like column order:

```text
ID1 hap1 ID2 hap2 Chr Start End LOD cM
```

See [`DATA_FORMAT.md`](DATA_FORMAT.md) for details.

## Metrics

Run:

```bash
python3 Metricas.py
```

Generated outputs:

```text
IBD_metrics_global.csv
IBD_metrics_per_pair.csv
IBD_metrics_by_length.csv
IBD_metrics_recall_accuracy_f1_by_length.pdf
```

In this analysis, `accuracy` is intentionally defined as:

```text
accuracy = TP / (TP + FP)
```

That is, it measures the fraction of LC-positive sequence that overlaps HC-positive sequence. This is equivalent to what is often called precision, but the output label is `accuracy` for this project.

Other definitions:

```text
TP = HC intersect LC
FN = HC intersect complement(LC)
FP = complement(HC) intersect LC
TN = complement(HC) intersect complement(LC)
recall = TP / (TP + FN)
f1 = 2TP / (2TP + FP + FN)
```

The length plot groups segments into `0.5 cM` bins and uses the bin midpoint on the x-axis.

## Individual IBD mosaic

Run:

```bash
Rscript MosaicoIBD.R
```

The target individual is set near the top of the script:

```r
base_ind <- "Bep03"
```

Generated outputs:

```text
IBD_mosaic_by_partner_Bep03.tiff
IBD_mosaic_by_partner_Bep03.csv
```

The mosaic shows one row per IBD partner and dataset. Segments are colored by length:

```text
Short  < 3 cM
Medium 3-5 cM
Long   >= 5 cM
```

## Notes

Large generated outputs and raw data directories are ignored by `.gitignore`. If you want to publish selected example figures, add them explicitly with `git add -f`.

## Circular summary plot

Run:

```bash
Rscript CircusIBD.R
```

The target pair is set near the top of the script:

```r
base_ind <- "Bep03"
partner_filter <- "LR174"
min_cm <- 0
```

Generated outputs:

```text
IBD_circus_Bep03_vs_LR174.tiff
IBD_circus_Bep03_vs_LR174.csv
```

The circular plot is intended as a compact pairwise summary. It shows two tracks, `hcWGS` and `lcWGS`, and colors segments by length class.
