# Identifiability of Contagion Components amid Environmental Fluctuations in Aggregated Default Counts

This repository contains the scientific code and workflows required to reproduce the numerical experiments and analysis for the manuscript "Identifiability of Contagion Components amid Environmental Fluctuations in Aggregated Default Counts".

## Overview

The repository provides the numerical engine to analyze aggregated annual default counts using hierarchical models. The project structure supports:
- Generating and processing aggregated annual default counts
- Evaluating static benchmark models (M0 NULL, single-factor Vasicek)
- Analyzing Probit-Normal environmental fluctuations (M1 / M2)
- Modeling contagion components with Gaussian-mixture models (M3)
- Cross-language validation (Python / Julia) to verify the cross-language validated best-found candidate
- Subperiod analysis

## Repository Structure

The core analysis is implemented in Jupyter Notebooks:

- `01_generate_synthetic_data.ipynb`: Generates a synthetic dataset for demonstration purposes.
- `02_numerical_experiments.ipynb`: Independent abstract numerical experiments.
- `03A_static_one_period_models.ipynb`: Independent static benchmark models.
- `03B_environmental_models.ipynb`: Independent environmental-model analysis; produces intermediate outputs for downstream M3 workflows.
- `03C_M3_python.ipynb`: Python M3 analysis.
- `03D_M3_julia.ipynb`: Julia M3 analysis.
- `03E_M3_cross_language_validation.ipynb`: Cross-language validation of Python and Julia M3 outputs.
- `03F_subperiod_analysis.ipynb`: Independent subperiod analysis.

### Directories and Modules

- `modules/`: Contains the canonical scientific Python implementation.
- `models.py`, `analysis_helpers.py`, `revision_models.py` (root level): Maintained alongside `modules/` to ensure exact backward compatibility and preserve canonical behavior for existing notebooks.
- `julia/`: Contains the Julia scientific implementation, dependencies (`Project.toml`, `Manifest.toml`), and test suites.
- `data/`: Contains the fully synthetic dataset for demonstration.
- `figure/`: Contains manuscript and diagnostic figures.
- `pdata/`: Local directory for generated intermediate outputs (not tracked in Git).
- `docs/`: Audits and implementation reports.

## Execution Workflows

### A. Exact Manuscript Reproduction

Exact reproduction of the manuscript numerical estimates requires the original proprietary dataset (`data/M_1920_2023.csv`). 

Execution order for manuscript reproduction:
1. `02_numerical_experiments.ipynb` (independent)
2. `03A_static_one_period_models.ipynb` (independent)
3. `03B_environmental_models.ipynb` (independent, produces local `pdata` outputs used downstream)
4. `03C_M3_python.ipynb` (depends on relevant `03B` outputs)
5. `03D_M3_julia.ipynb` (independent Julia M3 manuscript run using real data)
6. `03E_M3_cross_language_validation.ipynb` (cross-language validation requiring manuscript outputs/results from `03B`, `03C`, and `03D`)
7. `03F_subperiod_analysis.ipynb` (independent subperiod analysis reading real data directly; does NOT depend on `03B` or `pdata`)

*(Note: `01_generate_synthetic_data.ipynb` is not required for exact manuscript reproduction).*

### B. Public Synthetic Demonstration

The repository includes a synthetic demonstration workflow for users without access to the proprietary dataset.

**Warning: Synthetic results do NOT reproduce the manuscript numerical estimates.**

Execution order for synthetic demonstration:
1. `01_generate_synthetic_data.ipynb` (generates `data/M_1920_2023_synthetic.csv`)
2. `02_numerical_experiments.ipynb`
3. `03A_static_one_period_models.ipynb`
4. `03B_environmental_models.ipynb`
5. `03C_M3_python.ipynb`

*(Note: `03D` and `03E` strictly require the proprietary dataset and are excluded from the synthetic workflow. Additionally, while `03F_subperiod_analysis.ipynb` contains a synthetic data loader fallback, its full audit contains real-data-specific boundary checks. Therefore, exact manuscript reproduction of `03F` uses the proprietary dataset, and it is not claimed as a validated public synthetic demonstration.)*

## Data Availability / Privacy

Public repository includes:
- `data/M_1920_2023_synthetic.csv`

Public repository does NOT include:
- `data/M_1920_2023.csv`

The proprietary data cannot be redistributed. The public repository contains synthetic data only. The synthetic dataset is provided to demonstrate the code and workflow and does not reproduce the numerical estimates reported in the paper. Exact reproduction of the manuscript estimates requires the original same-schema proprietary dataset, which is not included in this repository.

## Intermediate Data (`pdata`) Policy

The `pdata/` directory acts as a locally generated intermediate and output directory and is explicitly excluded from the Git repository. During the manuscript workflow (e.g., `03B`, `03C`, `03E`), local `pdata` outputs may be produced and utilized by downstream notebooks. Proprietary-data-derived outputs stored in `pdata/` are never published.

## Figure Policy

The `figure/` directory contains manuscript figures and diagnostic plots. These are tracked by Git and available in the public release.

## Python Environment

The Python environment is managed via Conda. To set up the environment:

```bash
conda env create -f environment.yml
conda activate idms-jpsj
```

For batch execution or validation, you can invoke the canonical environment directly:

```bash
conda run -n idms-jpsj python [script.py]
conda run -n idms-jpsj jupyter nbconvert --ExecutePreprocessor.kernel_name=idms-jpsj [notebook.ipynb]
```

## Julia Environment

The Julia component utilizes Julia 1.12.7 as the canonical version. The environment configuration is defined in:
- `julia/Project.toml`
- `julia/Manifest.toml`

To set up the Julia environment:

```bash
cd julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The Julia implementation includes a canonical test target where 57 / 57 tests PASS. For more details on running optimizations and tests, please see `julia/README.md`.

## License and Citation

This software is released under the [MIT License](LICENSE).

If you use this software, please cite the associated paper:

Shintaro Mori,
"Identifiability of Contagion Components amid Environmental Fluctuations in Aggregated Default Counts,"
arXiv:2604.18118, 2026.
https://doi.org/10.48550/arXiv.2604.18118

BibTeX:
```bibtex
@article{mori2026identifiability,
  author = {Shintaro Mori},
  title = {Identifiability of Contagion Components amid Environmental Fluctuations in Aggregated Default Counts},
  journal = {arXiv preprint arXiv:2604.18118},
  year = {2026},
  doi = {10.48550/arXiv.2604.18118}
}
```
