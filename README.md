# Contagion vs Macroeconomic Fluctuations in Default Data

Code and data for the paper:

**Shintaro Mori,  
[*Contagion or Macroeconomic Fluctuations? Identifiability in Aggregated Default Data*](https://arxiv.org/abs/2604.18118)**  
arXiv:2604.18118

## Overview

This repository contains the code used to analyze default clustering
and to decompose aggregate fluctuations into contagion and macroeconomic components.

The main objective is to assess whether contagion effects can be
identified separately from macroeconomic fluctuations using aggregated
annual default count data.

The repository is intended to ensure reproducibility of the results reported in the paper.

## Paper

- **arXiv**: 2604.18118
- **Title**: *Contagion or Macroeconomic Fluctuations? Identifiability in Aggregated Default Data*
- **Author**: Shintaro Mori

## Contents

- `IDM_Data_Analysis.ipynb`  
  Empirical analysis of default data, including likelihood estimation,
  variance decomposition, and tail-risk comparison.

- `IDM_Simulation.ipynb`  
  Simulation of contagion models and validation of decomposition results.

## Reproducibility

The main results of the paper can be reproduced as follows:

1. Run `IDM_Data_Analysis.ipynb`  
   → reproduces parameter estimation, likelihood comparison, variance decomposition,  
   and the empirical figures and tables.

2. Run `IDM_Simulation.ipynb`  
   → reproduces simulation-based validation, distributional comparisons,  
   and robustness checks.

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
- The repository is designed to accompany the arXiv paper and support reproducibility.

## Data and code availability

This repository contains the scripts used to generate the figures
and tables reported in the paper.

## Citation

If you use this code, please cite:

```bibtex
@article{mori2026contagion,
  author  = {Shintaro Mori},
  title   = {Contagion or Macroeconomic Fluctuations? Identifiability in Aggregated Default Data},
  journal = {arXiv preprint arXiv:2604.18118},
  year    = {2026}
}
