# FRET_TEVp_26
This repository contains a pipeline used to process fluorometric data produced by the FRET assay described in section 2.2 in "Impact of P1' substrate substitutions on the catalytic efficiency of Tobacco Etch Virus Protease"

It is recommended to run the R and Rmd files in R studio, and the Jupyter notebook activity_analysis.ipynb in VS Code or another Jupyter notebook interpreter such as Google Colab.

The original fluorometric data can be found in the file: 'SUBGVP_Esin_20feb26.xlsx'

The plate map for the assay can be found in the file: 'SUBGVP_Esin_20feb26_platemap.csv'

## Acknowledgements
The pipeline is adapted from Marius Krogsgaard Thomsen's original pipeline of the FRET-based kinetic assay. This pipeline is adapted to accomodate multiple substrates, whilst the original pipeline was constructed to analyse multiple TEVp variants.

## Description of data pipeline, including modules and packages
The pipeline contains the following files:
- parse_plate.R: to clean the original fluorometric data
- core_analysis_v2.R: to correct for photobleaching, normalize the data and calculate initial rates
- activity_analysis.ipynb: to find the $k_obs$ and plotting.

To run parse_plate.R and core_analysis_v2.R, use the driver FRET_Driver.Rmd.
The R files contains various which can be downloaded in the first code block of FRET_Driver.RmD:
```{r}
install.packages('optparse')
install.packages('readxl')
install.packages('dplyr')
install.packages('tidyr')
install.packages('stringr')
install.packages('lubridate')
install.packages('purrr')
install.packages('readr')
install.packages('writexl')
install.packages('stats')
```

The file activity_analysis.ipynb uses the following packages:
- pandas 
- matplotlib
- numpy
- scipy
These can be downloaded in the terminal as follows:
```{bash}
pip install [insert package name]
```

## How to run the pipeline
0. Make sure to have an excel sheet of fluorometric data and csv file with plate map
1. Run FRET_Driver.Rmd in order to run plate_map.R and core_analysis_v2.R
2. Take the output QC file and run it in activity_analysis.ipynb
    - Due to the many plots and tables produced in activity_analysis.ipynb, it is recommended to separate folders for them.

An example of how to run the pipeline can be found in the folder 'Output'
