# Stage 4A Audit: Notebook 03B Dependency and Import Safety

## 1. Environment Preflight
- **Canonical conda environment**: `idms-jpsj`
- **Python executable**: `[LOCAL_PATH]`
- **Python version**: 3.11.16
- **Jupyter kernel**: `idms-jpsj`

## 2. Git Preflight
- **Branch**: `code-publication-refactor`
- **Status**: Working tree clean
- **Commit**: `e6edd70 Migrate static benchmark notebook to shared modules`

## 3. All Function Inventory and 4. Exact/Different Comparison
| Function Name | Location | Module | Status | AST Hash |
|---|---|---|---|---|
| log_comb | Cell 3 | models | EXACT | f5549af3 |
| logsumexp_arr | Cell 3 | models | EXACT | 71b2307f |
| logsumexp | Cell 3 | N/A | NOT IN MODULES | a455990a |
| logsubexp | Cell 3 | models | EXACT | 34c0f504 |
| normalize_logpmf | Cell 3 | models | EXACT | 46031ac4 |
| var_from_pmf | Cell 3 | models | EXACT | 913f8734 |
| es_from_pmf | Cell 3 | models | EXACT | 1382bcbd |
| moments_from_pmf | Cell 3 | models | EXACT | 90e096dc |
| m_rho_from_pmf | Cell 3 | models | EXACT | fce6b375 |
| survival_from_pmf | Cell 3 | models | EXACT | 8797edaa |
| tail_metrics_from_pmf | Cell 3 | models | EXACT | 9d8205d5 |
| pmf_Torri | Cell 3 | models | EXACT | 383567ce |
| moments_from_params_Torri | Cell 3 | models | DIFFERENT | ab3a9aa5 |
| pmf_LD | Cell 3 | models | EXACT | 94bdebcd |
| m_rho_LD | Cell 3 | models | EXACT | f59121ec |
| check_mrho_formula_vs_pmf_LD | Cell 3 | models | EXACT | 89748549 |
| phi2_aa | Cell 3 | models | EXACT | dce57d04 |
| default_corr_from_asset_corr | Cell 3 | models | EXACT | 8a8958c8 |
| invert_asset_corr_from_default_corr | Cell 3 | models | EXACT | 4f21e66a |
| pmf_Vas | Cell 3 | models | EXACT | 7d395b89 |
| log_pmf_Torri | Cell 4 | models | EXACT | 85119dfa |
| logsumexp_vec | Cell 5 | models | EXACT | b7e5b669 |
| log_pmf_LD_1pt | Cell 5 | models | EXACT | a5b02e29 |
| sigmoid | Cell 7 | analysis_helpers | EXACT | 0c8583d8 |
| softplus | Cell 7 | analysis_helpers | EXACT | 02e67d46 |
| inv_softplus | Cell 7 | analysis_helpers | EXACT | e19a2503 |
| _clean_nh | Cell 7 | analysis_helpers | EXACT | 85ca9f9e |
| make_grouped_nh | Cell 7 | analysis_helpers | EXACT | 6e2f3833 |
| _clean_nL | Cell 7 | analysis_helpers | EXACT | f3850b0c |
| scaled_counts | Cell 7 | analysis_helpers | EXACT | c849991f |
| gh_normal_weights | Cell 7 | analysis_helpers | EXACT | b452d24b |
| aggregate_m_q_to_rho | Cell 7 | analysis_helpers | EXACT | 4c2102ac |
| torri_m_q_for_n | Cell 7 | analysis_helpers | EXACT | 40b3357e |
| build_compare_decomp_table | Cell 7 | analysis_helpers | DIFFERENT | 3c84b521 |
| ld_m_q_for_n | Cell 7 | analysis_helpers | EXACT | e965c611 |
| log_pmf_pn_LD_1pt | Cell 9 | N/A | NOT IN MODULES | da055f67 |
| neg_ll_pn_LD | Cell 9 | N/A | NOT IN MODULES | 4375a0c4 |
| fit_one_pn_LD | Cell 9 | N/A | NOT IN MODULES | 5b45dbeb |
| posterior_mean_pt_pn_LD | Cell 9 | N/A | NOT IN MODULES | 114a778f |
| add_posterior_mean_pt_column_pn_LD | Cell 9 | N/A | NOT IN MODULES | c502fbbb |
| log_pmf_pn_Torri_1pt | Cell 11 | N/A | NOT IN MODULES | 2316c073 |
| neg_ll_pn_Torri | Cell 11 | N/A | NOT IN MODULES | 1d94244f |
| neg_ll_pn_null | Cell 11 | N/A | NOT IN MODULES | bdfd66fa |
| fit_pn_null_from_grouped | Cell 11 | N/A | NOT IN MODULES | ba92bd9b |
| fit_one_pn_Torri | Cell 11 | N/A | NOT IN MODULES | 5da3c4da |
| posterior_mean_pt_pn_Torri | Cell 11 | N/A | NOT IN MODULES | 633f3fbe |
| add_posterior_mean_pt_column_pn_Torri | Cell 11 | N/A | NOT IN MODULES | fb798acc |
| print_pn_torri_fit | Cell 11 | N/A | NOT IN MODULES | 845c08bd |
| variance_decomp_pn_ld_fixed_nbar | Cell 17 | N/A | NOT IN MODULES | 6645b4e8 |
| variance_decomp_pn_torri_fixed_nbar | Cell 18 | N/A | NOT IN MODULES | 43156f55 |
| empirical_pmf_from_counts | Cell 22 | analysis_helpers | DIFFERENT | e0ec46fd |
| pmf_pn_null_fixed_n | Cell 23 | N/A | NOT IN MODULES | 9a5c34d3 |
| pmf_pn_ld_fixed_n | Cell 23 | N/A | NOT IN MODULES | 8a524006 |
| pmf_pn_torri_fixed_n | Cell 23 | N/A | NOT IN MODULES | 2eb0910a |
| plot_survival_m012_one | Cell 24 | N/A | NOT IN MODULES | ffeb388f |
| plot_pmf_m012_one | Cell 25 | N/A | NOT IN MODULES | fdaa1ca3 |
| pmf_ld_sparse_full | Cell 27 | N/A | NOT IN MODULES | 45297537 |
| pmf_m1_full | Cell 27 | N/A | NOT IN MODULES | a0f03584 |
| log_pmf_torri_full_vec | Cell 27 | N/A | NOT IN MODULES | abac5e68 |
| pmf_m2_full | Cell 27 | N/A | NOT IN MODULES | 8f15be23 |
| pmf_m0_full | Cell 27 | N/A | NOT IN MODULES | 52341872 |
| directed_kl | Cell 27 | N/A | NOT IN MODULES | f92c40ca |
| mean_kl_for_class | Cell 27 | N/A | NOT IN MODULES | f3182392 |

## 5. Global/Free-Name Analysis & 6. Import-safe & 7. Retain-inline
Functions like `log_pmf_Torri` (depends on DIFFERENT `logsumexp`), `log_pmf_LD_1pt` (depends on `logsumexp_vec`), and their nested callers (`neg_ll_pn_Torri`, `fit_one_pn_Torri`, etc.) use notebook-specific overrides.
- **Import-safe candidates (no conflicting global deps):**
  - `inv_softplus`, `softplus`, `sigmoid`, `_clean_nh`, `_clean_nL`, `gh_normal_weights`, `scaled_counts`, `torri_m_q_for_n`, `ld_m_q_for_n`, `var_from_pmf`, `es_from_pmf`, `moments_from_pmf`, `m_rho_from_pmf`, `pmf_Torri`, `pmf_Vas`, `default_corr_from_asset_corr`, `invert_asset_corr_from_default_corr`, `check_mrho_formula_vs_pmf_LD`, `tail_metrics_from_pmf`
- **Retain-inline candidates:**
  - `logsumexp`, `log_pmf_Torri`, `log_pmf_LD_1pt`, `moments_from_params_Torri`, `m_rho_LD`, `build_compare_decomp_table`, `empirical_pmf_from_counts`, `log_comb` and all other `NOT IN MODULES` fitting wrappers.

## 8. 03B / 03F / revision_models Comparison
The 10 exact duplicate functions between 03B and 03F (`fit_one_pn_LD`, `neg_ll_pn_Torri`, etc.) **DO NOT EXIST** by name in `revision_models.py`. The `revision_models.py` uses different names/implementations (e.g., `nll_probit_normal`, `fit_one_hier_Torri_multistart`). They must remain inline.

## 9. `make_grouped_nh` Comparison
The 03B version of `make_grouped_nh`:
- Groups using `tmp.value_counts().reset_index(name="cnt")` (which sorts by frequency descending, not mathematically sorting by `n` and `h`).
- Dtypes are strictly `int`.
- Excludes `n <= 0`, `h < 0`, and `h > n` internally via `_clean_nh`.
Since grouping convention impacts fitting (order of operations in log-sums), this must remain inline or uniquely tracked.

## 10. M0/M1/M2 Dependency Graph
- **M0 (Probit-Normal Null)**: `fit_pn_null_from_grouped` -> `neg_ll_pn_null` -> `gh_normal_weights`, `norm.cdf`, `binom.logpmf`
- **M1/M2 Torri**: `fit_one_pn_Torri` -> `make_grouped_nh`, `neg_ll_pn_Torri` -> `log_pmf_pn_Torri_1pt` -> `log_pmf_Torri` -> `logsumexp`, `logsubexp`, `log_comb`
- **M1/M2 LD**: `fit_one_pn_LD` -> `make_grouped_nh`, `neg_ll_pn_LD` -> `log_pmf_pn_LD_1pt` -> `logsumexp_vec`, `log_comb`

## 11. Canonical-result-sensitive Functions
The functions governing the M0 null hypothesis and M2 boundaries (`v=0` logic):
- `fit_pn_null_from_grouped`
- `fit_one_pn_Torri` (contains the `v_raw < v_zero_tol` fallback logic)
- `neg_ll_pn_Torri`
- `log_pmf_Torri` (and its nested dependencies like `logsumexp`)

## 12. Namespace Collision Audit
- The following names are used as variables in a `Store` context within the notebook: `m`, `model`, `models`, `mu`, `p`, `q`, `sigma`.
- **Conclusion**: We CANNOT use `import modules.models as m`, `import modules.models`, or `import modules.models as models`. We MUST use explicit `from modules.x import y`.

## 13. Output/Write Inventory
- **Cell 15**: `to_csv`
- **Cell 19**: `to_csv`
- **Cell 20**: `savefig`
- **Cell 23**: `to_csv`
- **Cell 24**: `savefig`
- **Cell 25**: `savefig`
- **Cell 27**: `to_csv`
- Generates `pdata` files which are consumed by plotting cells or downstream notebooks.

## 14. Data Mode Audit
- Notebook explicitly checks `if REAL_DATA.exists():` and explicitly alerts the user with `DATA SOURCE: FULLY SYNTHETIC ... The resulting numerical estimates will NOT reproduce the manuscript` if it falls back. The synthetic fallback is **EXPLICIT**, not silent.

## 15. Recommended Stage 4B Minimal Migration Set
Only truly independent, import-safe functions with zero dependency on notebook-redefined globals:
`sigmoid`, `softplus`, `inv_softplus`, `_clean_nh`, `_clean_nL`, `gh_normal_weights`, `scaled_counts`, `torri_m_q_for_n`, `ld_m_q_for_n`, `var_from_pmf`, `es_from_pmf`, `moments_from_pmf`, `m_rho_from_pmf`, `pmf_Torri`, `pmf_Vas`, `default_corr_from_asset_corr`, `invert_asset_corr_from_default_corr`, `check_mrho_formula_vs_pmf_LD`, `tail_metrics_from_pmf`.

## 16. Functions that MUST NOT be migrated
- `logsumexp`, `moments_from_params_Torri`, `m_rho_LD`, `log_pmf_Torri`, `log_pmf_LD_1pt`
- All fitting wrappers (`fit_one_pn_Torri`, `fit_one_pn_LD`, `neg_ll_pn_Torri`, `make_grouped_nh`, etc.)
