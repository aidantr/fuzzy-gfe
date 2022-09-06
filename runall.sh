#!/bin/bash

# Runall for "A Fuzzy Clustering Approach to Estimating Grouped Fixed-Effects" by Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)
# Runs Matlab files to produce all tables and figures in the main text and appendix

# set main directory
cd <put path here>

# run Matlab files in order
matlab -batch "run('code/fig1.m')"
matlab -batch "run('code/fig2.m')"
matlab -batch "run('code/simulate_panel.m')"
matlab -batch "run('code/table1.m')"
matlab -batch "run('code/fig3.m')"
matlab -batch "run('code/tableB1.m')"


