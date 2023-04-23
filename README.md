# SMG6 RNA-Seq Data Analysis

This repository contains the code and results of an analysis of a bulk RNA-Seq dataset originally published by Guerra et al. (on Nov 29, 2021), which can be downloaded from [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8699217/). 

The raw data was first processed using shell scripts, then analyzed using R. The details can be found in the "report_and_tutorial" folder. The folder contains a detailed tutorial named "tutorial.Rmd" on how to run the scripts to process the raw data and analyze the sequencing results. The tutorial file can be opened in RStudio.

## Folder Structure
* report_and_tutorial: This folder contains the tutorial file (tutorial.Rmd) that explains the analysis steps in detail and a project report (report.Rmd). 
* scripts: This folder contains the shell scripts used for processing the raw data.
* outputs: This folder contains the multiqc outputs.
* key_datasets: the important datasets generated from the analysis. 
* images: images used in the report and tutorial. 

## Conda environment 

```bash
conda env create --name myenv --file angsd_project_environment.yml
```