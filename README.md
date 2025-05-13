Replication Package for
# “Selective Emigration and Economic Development: A Generalized Approach”

_Narcisse Cha’ngom, Christoph Deuster, Frédéric Docquier,_ and _Joël Machado_


Version: May 2025


# Overview
This replication package includes all data and code required to reproduce the tables and figures presented in the paper. The datasets used in the paper are provided. The replication files are organized into four folders and three master scripts:

• appendix

• main

• matlab

• stata

• runme


# Stata
This folder contains a single subfolder named “data”, which includes the cleaned datasets used for the empirical analysis.

#  Matlab
This folder contains three subfolders:
  1. data – Includes the input data for model simulations, such as the bilateral migration matrix and country-level characteristics.
  2. scripts – Contains all MATLAB (.m) files used to solve the theoretical model and perform all the counterfactual simulations.
  3. output – Stores the output files from the simulations, including those presented in the paper and online appendix.

# Main
This folder includes:
  1. A subfolder 'figures' and another subfolder 'tables', both replicating all figures and tables in the main text.
  2. A Stata do-file named _main_results.do_, which replicates all main text results.

# Appendix
This folder includes:
  1. Subfolders 'figures' and 'tables' containing all figures and tables in the online appendix.
  2. A Stata do-file named _appendix_results.do_, which replicates these appendix materials.

# Runme
These are master scripts that execute the full set of simulations and empirical analyses. The should be executed in the following order:
- _runme_main.m_: MATLAB script that runs all model simulations and exports the benchmark results.
- _runme_sensitivity.m_: MATLAB script that tests sensitivity of the benchmark simulations to alternative parameter values.
- _runme.do_: Stata master do-file for both empirical estimations and simulation-related procedures.
Software
The empirical analysis was conducted using Stata MP Version 18. The theoretical model simulations were implemented in Matlab Version 2024b.
