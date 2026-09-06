# Stage 2 Regression Audit: Notebook 02

- **baseline tag:** canonical-pre-refactor-20260906
- **refactor branch:** code-publication-refactor
- **canonical environment:** idms-jpsj
- **Python executable:** [LOCAL_PATH]
- **Python version:** 3.11.16
- **Jupyter kernel:** idms-jpsj

## Execution Method
Baseline and current versions of `02_numerical_experiments.ipynb` were both executed top-to-bottom in the identical canonical environment to ensure consistency.

## Numerical Comparison Results
- **numerical/text/table outputs:** 0 mismatches

## Refactor Scope
- **Stage 2 imported EXACT-AST functions:**
  - `default_corr_from_asset_corr`
  - `es_from_pmf`
  - `find_params_on_iso_Torri`
  - `invert_asset_corr_from_default_corr`
  - `m_rho_from_pmf`
  - `moments_from_pmf`
  - `pmf_Torri`
  - `pmf_Vas`
  - `tail_metrics_from_pmf`
  - `var_from_pmf`

## Notes
- `b_of_s` and `bisect_root` are nested helpers within `find_params_on_iso_Torri` and are not standalone import targets.
- The local variable `m` was found to be used as an existing variable inside the notebook. Therefore, `import modules.models as m` was not used. Instead, explicit `from modules.models import (...)` was employed.
- No DIFFERENT functions were modified.
- Figure binary files were not committed because plotting logic remains unchanged. Only the underlying numerical data was confirmed to be identical.

## Conclusion
- **Stage 2 conclusion:** PASS
