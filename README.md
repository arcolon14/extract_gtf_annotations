# extract_gtf_annotations
Extract target annotation from a GTF based on a set of target (outlier) sites.

## Usage

```sh
outliers_to_annotations.py started on 2024-01-31 21:30:08

usage: outliers_to_annotations.py [-h] -g GTF -l OUTLIER_LIST [-o OUTDIR] [-b BUFFER]

Subset an annotation GTF file to extract annotations with a specified distance from a set of target
outlier SNPs.

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     (str) Path to annotation GTF
  -l OUTLIER_LIST, --outlier-list OUTLIER_LIST
                        (str) Path to outlier list CSV
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory
  -b BUFFER, --buffer BUFFER
                        (int) Distance in bp plus/minus target SNP to extract a gene
                        [default=10000]
```
