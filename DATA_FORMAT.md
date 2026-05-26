# Data format

The scripts expect one gzipped IBD file per chromosome and dataset.

## Directory layout

```text
HC/chr1_ibd.ibd.gz
HC/chr2_ibd.ibd.gz
...
LC/chr1_ibd.ibd.gz
LC/chr2_ibd.ibd.gz
...
```

## IBD columns

Each row must contain at least nine whitespace-separated columns:

```text
ID1 hap1 ID2 hap2 Chr Start End LOD cM
```

| Column | Description |
|---|---|
| `ID1` | First sample ID. Coverage suffixes such as `_1X` and `_30X` are removed during metric comparison. |
| `hap1` | Haplotype label for `ID1`. This is retained in the raw data but not used as a matching key in the final metric script. |
| `ID2` | Second sample ID. |
| `hap2` | Haplotype label for `ID2`. |
| `Chr` | Chromosome name, expected as `chr1` to `chr22`. |
| `Start` | Segment start coordinate in base pairs. |
| `End` | Segment end coordinate in base pairs. |
| `LOD` | RefinedIBD LOD score. |
| `cM` | Segment length in centimorgans. |

Example:

```text
Bep03_30X 1 SAL111_30X 1 chr1 4294402 4955407 27.31 1.745
```

## Metric matching key

For the final metric script, segments are matched by:

```text
chromosome + sample1 + sample2
```

Sample order is normalized, so `A-B` and `B-A` are treated as the same pair. Haplotype labels are not part of the final matching key because phase switches can make haplotype-specific matching too strict.

## RefinedIBD outputs

The metric and mosaic scripts currently read the standard RefinedIBD files:

```text
chr<CHR>_ibd.ibd.gz
```

The SLURM scripts also produce merged files:

```text
chr<CHR>_ibd_merged.ibd
```

Those merged files are useful for downstream analyses only if they contain actual segment rows. If the file starts with `usage:`, the merge command was run incorrectly and the file should be regenerated.
