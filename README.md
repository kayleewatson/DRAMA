# DRAMA
Differential RNA Modification Analysis (DRAMA) using Oxford Nanopore direct RNA sequencing

Required packages:

* ggplot2
* ggrepel
* zoo
* ggpubr
* tidyverse
* dplyr
* tidyr
* optparse

```
Usage: DRAMA_dev_both_strands_density.r [options]


Options:
        -x REFSTATS1, --refstats1=REFSTATS1
                Refstats (mean ionic current) file for sample 1 [required]

        -y REFSTATS2, --refstats2=REFSTATS2
                Refstats (mean ionic current) file for sample 2 [required]

        -k KSINPUT, --ksinput=KSINPUT
                Refstats KS test file [required]

        -p PVALUE, --pvalue=PVALUE
                Cutoff for generalized ESD p-value [default 0.05]

        -c KSNORM, --ksnorm=KSNORM
                Cutoff value for the normalized KS statistic [default 3]

        -n NAME, --name=NAME
                Chromosome name or file prefix [default DRAMA]

        -l, --labels
                Include top 10 outlier labels on the plot (default is FALSE)

        -d HEIGHT, --height=HEIGHT
                User-provided height for each plot (default=5)

        -w WIDTH, --width=WIDTH
                User-provided width for each plot (default=6)

        -h, --help
                Show this help message and exit
```
