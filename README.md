# Replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)

This repository contains code to replicate the results of  "Approximating Grouped Fixed Effects Estimation via Fuzzy Clustering Regression." The paper is available [here](https://drive.google.com/file/d/1U_MJHtJcB7H1Edv3xceilU_HJoxhLssP/view). All files are written in Matlab.

## Required Matlab toolboxes

- Statistics and machine learning toolbox
- Optimization toolbox

## How to run

The bash file `runall.sh` runs all the necessary Matlab files to reproduce the tables and figures in the main text and appendix. If you're interested in using the FCR estimator for your own dataset, check out the `FCR.m` function in the code/functions folder.

## Repository structure

- `data/raw` contains the raw data, which come from the Bonhomme and Manresa (2015) [replication files](https://www.dropbox.com/s/ssjabvc2hxa5791/Bonhomme_Manresa_codes.zip?dl=0)
- `data/intermediate` stores intermediate files
- `code` contains the Matlab functions needed to implement our estimator and all files to produce the output
- `output` stores results for all tables and figures in the paper

## Parallelization

Our main estimation is run with 250 parallel cores. However, the code can be run with any number of cores (just adjust `parpool`) although this will change computation time.






