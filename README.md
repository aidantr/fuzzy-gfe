# Replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)

This repository contains all the code to replicate the results of  "A Fuzzy Clustering Approach to Estimating Grouped Fixed-Effects." The paper is available [put link here]. All the files are written in Matlab.

## Required Matlab toolboxes

- Statistics and machine learning toolbox
- Optimization toolbox
- Deep learning toolbox

## How to run

The bash file `runall.sh` runs all the necessary Matlab files in order to reproduce the tables and figures in the paper. 

## Repository structure

- `data/raw` contains the raw data, which come from the Bonhomme and Manresa (2015) replication files
- `data/intermediate` stores intermediate files
- `code` contains the Matlab functions needed to implement our estimator and all files to produce the output
- `output` stores results for all tables and figures in the main text and appendix

## Parallelization

Our main estimation is run with 250 parallel cores. However, the code can be run with any number of cores (just adjust `parpool(cores)`) although this will change computation time.






