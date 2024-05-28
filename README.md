# scRNA-seq/scATAC-seq working scripts

This repository holds the two main R scripts used in this study.

## scRNA-seq

This script is called _scRNAseq.Rmd_. It is a **R Markdown** script that outputs an interactive HTML file with all graphs and results. It performs the full **SCENIC** pipeline on the data, so as you can imagine, this is **EXTREMELY** time and memory-intensive.

On average, the analysis took around 1 and 2 days and the maximum memory usage was around 30 and 40 GB of RAM.

These are the main steps:
* Place the script in the same folder with the data to analyze.
* Modify the appropriate variables within it, namely the data files and such.
* Once you think you are ready, launch it with

```{bash}
Rscript -e "require(rmarkdown); rmarkdown::render(scRNAseq.Rmd, quiet=TRUE)"
```

Ensure to have RMarkdown installed (R):

```{R}
install.packages('rmarkdown')
```

Once the main SCENIC pipeline starts, if it ever fails or is interrupted in any of its 4 steps, it should be able to continue from that point if the script is executed again.

If all goes well, at the end you should have a HTML file of around 15 MB in your hands with all your results. Executing the script again will **SKIP** all the SCENIC pipeline because it will know it successfully finished. The file responsible for this behaviour is placed usually in _int/scenicOptions.Rds_. If you don't ever modify it, it lets you execute the script again and waste no time repeating things. This lets you, for example, adjust some graphical analysis done after SCENIC without having to call SCENIC again.

## scATAC-seq

This script is called _scATAC-seq.R_. It is a **R** file. This file is easier on the machine, and also easier to execute:

```{bash}
R -f scATAC-seq.R
```

Ensure that it is placed in the same directory as your data. Its output is in the form of a **PDF** file with all the graphical results. It is less memory and time-intensive, although I would expect no less than 20 GB or RAM to be necessary, depending on your data.

Some checkpoints were placed in order to check the parts that fail more frequently. This way, since these parts are lengthy, once you get them right they are skipped if the script is executed again.

## Notes

* The garbage collector is called several times in both scripts (`gc()`) in order to get rid of variables that are not used again
* The preferred format for saving files is `RDS`. Load them with `loadRDS()` and save them with `saveRDS()`
* The performance reported was tested in a Lenovo Ideapad 330 with an AMD Ryzen 3 with 4 cores, 4 threads, running Linux and with 12 GB of RAM, artificially expanded through swap space to 62 GB in an SSD.* 

## Utilities
* BioNERO: https://github.com/almeidasilvaf/BioNERO
* SCENIC: https://github.com/aertslab/SCENIC
* SCOPELOOMR: https://github.com/aertslab/SCopeLoomR
* AUCell: https://github.com/aertslab/AUCell
* GOFuncR: https://github.com/sgrote/GOfuncR
* SEURAT: https://github.com/satijalab/seurat
* SIGNAC: https://github.com/stuart-lab/signac
