# Stage 3 Pre-Implementation Audit: Notebook 03A

## Audit Results

- **Candidate**: summarize_nt_default
  - **Notebook location**: Cell 3
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 02f32a9e
  - **Status**: EXACT
  - **Free/global names**: L, df, float, label, n, np
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: summarize_nt_default_by_period
  - **Notebook location**: Cell 4
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 235623e2
  - **Status**: EXACT
  - **Free/global names**: L, df, float, g, label, n, np, out, period
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: fit_one_Torri
  - **Notebook location**: Cell 10
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 7d603fb2
  - **Status**: EXACT
  - **Free/global names**: bool, data, float, h, int, ll, log_pmf_Torri, lp, math, maxiter, minimize, n, neg_ll_torri, np, p, res, series_h, series_n, str, th, theta, u, v, x0
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: neg_ll_torri
  - **Status**: NOT FOUND as top-level function in Notebook. Likely a nested helper.

- **Candidate**: fit_one_LD
  - **Notebook location**: Cell 12
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 417a0e77
  - **Status**: EXACT
  - **Free/global names**: bool, c, c_arr, float, g, h, h_arr, int, ll, log_pmf_LD_1pt, lp, math, minimize, n, n_arr, neg_ll_ld, np, p, pd, q, res, series_h, series_n, str, th, theta, tmp, x0, zip
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: neg_ll_ld
  - **Status**: NOT FOUND as top-level function in Notebook. Likely a nested helper.

- **Candidate**: torri_m_q_fn
  - **Notebook location**: Cell 18
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 96fbc21e
  - **Status**: EXACT
  - **Free/global names**: _fn, float, m_n, moments_from_params_Torri, n, np, p, q_n, rho_n, u, v
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: ld_m_q_fn
  - **Notebook location**: Cell 18
  - **Module location**: modules.analysis_helpers
  - **AST hash**: c6936fd3
  - **Status**: EXACT
  - **Free/global names**: _fn, float, m_n, m_rho_LD, n, np, p, q, q_n, rho_n
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: vas_m_q_fn
  - **Notebook location**: Cell 18
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 8e949e87
  - **Status**: EXACT
  - **Free/global names**: _fn, a, float, norm, p, phi2_aa, q_pair, rhoA
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: to_float_dict
  - **Notebook location**: Cell 16
  - **Module location**: modules.analysis_helpers
  - **AST hash**: 8fc3d9bc
  - **Status**: EXACT
  - **Free/global names**: Exception, d, float, k, out, v
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: log_pmf_Torri
  - **Notebook location**: Cell 10
  - **Module location**: modules.models
  - **AST hash**: 85119dfa
  - **Status**: EXACT
  - **Free/global names**: a, ab, b, bc, c, diff_log, float, h, inner, int, la, log_comb, logsubexp, logsumexp, lp, math, n, np, p, pb, t0, t1, u, v
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

- **Candidate**: log_pmf_LD_1pt
  - **Notebook location**: Cell 12
  - **Module location**: modules.models
  - **AST hash**: a5b02e29
  - **Status**: EXACT
  - **Free/global names**: eps, float, h, k, log_1p, log_1q, log_1r, log_bin, log_comb, log_p, log_pk, log_r, logsumexp_vec, math, n, nk, np, p, q, range, rk, t, terms
  - **Global dependency equivalence**: OK (no suspicious globals)
  - **Import-safe**: YES

## Write Inventory
- **Writes found**: `savefig`, `to_csv`, `pdata access`

## Import Alias Check
- **Aliases used as variables**: `m`
- **Action**: Do not use `import modules.models as m`. Use explicit `from modules.models import (...)`.
