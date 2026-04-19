# Contagion vs Macroeconomic Fluctuations in Default Data

Code and data for the paper:

**"Contagion or Macroeconomic Fluctuations? Identifiability in Aggregated Default Data"**

## Overview

This repository contains the code used to analyze default clustering
and to decompose aggregate fluctuations into contagion and macroeconomic components.

The main objective is to assess whether contagion effects can be
identified separately from macroeconomic fluctuations using aggregated
annual default count data.

The repository is intended to ensure reproducibility of the results reported in the paper.

## Contents

- `IDM_Data_Analysis.ipynb`  
  Empirical analysis of default data, including likelihood estimation
  and variance decomposition.

- `IDM_Simulation.ipynb`  
  Simulation of contagion models and validation of decomposition results.

## Reproducibility

The main results of the paper can be reproduced as follows:

1. Run `IDM_Data_Analysis.ipynb`  
   → reproduces parameter estimation, likelihood comparison, and variance decomposition.

2. Run `IDM_Simulation.ipynb`  
   → reproduces simulation-based validation and robustness checks.

All figures and tables in the paper are generated from these scripts.

## Requirements

- Python 3.x
- Standard scientific libraries:
  - numpy
  - scipy
  - pandas
  - matplotlib

## Notes

- The analysis is based on aggregated annual default count data.
- The focus is on structural interpretation rather than predictive performance.

## Data and code availability

This repository contains all scripts used to generate the figures
and tables in the paper.

## Reference

If you use this code, please cite:

Mori, S.  
*Contagion or Macroeconomic Fluctuations? Identifiability in Aggregated Default Data.*

## License

This project is licensed under the MIT License.
