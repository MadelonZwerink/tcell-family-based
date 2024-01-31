# T cell family based

This project shows the sizes of different T cell families.

## Usage

Download the folder and run the run_model file.

## Project Structure

The project structure distinguishes three kinds of folders: - read-only (RO): not edited by either code or researcher - human-writeable (HW): edited by the researcher only. - project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.

```         
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── requirements.txt
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── R                  <- Source code for this project (HW)
```

## Session info

```         
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=Dutch_Netherlands.utf8  LC_CTYPE=Dutch_Netherlands.utf8   
[3] LC_MONETARY=Dutch_Netherlands.utf8 LC_NUMERIC=C                      
[5] LC_TIME=Dutch_Netherlands.utf8    

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] plyr_1.8.9        cowplot_1.1.3     data.table_1.15.0 magrittr_2.0.3   
 [5] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4      
 [9] purrr_1.0.2       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1     
[13] ggplot2_3.4.4     tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] gtable_0.3.4      compiler_4.3.2    renv_1.0.3        tidyselect_1.2.0 
 [5] Rcpp_1.0.11       scales_1.3.0      R6_2.5.1          labeling_0.4.3   
 [9] generics_0.1.3    munsell_0.5.0     pillar_1.9.0      tzdb_0.4.0       
[13] rlang_1.1.1       utf8_1.2.4        stringi_1.7.12    timechange_0.3.0 
[17] cli_3.6.1         withr_3.0.0       grid_4.3.2        rstudioapi_0.15.0
[21] hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       glue_1.6.2       
[25] farver_2.1.1      fansi_1.0.6       colorspace_2.1-0  tools_4.3.2      
[29] pkgconfig_2.0.3  
```

## Add a citation file

Create a citation file for your repository using [cffinit](https://citation-file-format.github.io/cff-initializer-javascript/#/)

## License

This project is licensed under the terms of the [MIT License](/LICENSE).
