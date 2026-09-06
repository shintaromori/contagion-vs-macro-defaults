# Code Publication Refactor Implementation Plan

Created: 2026-09-06  
Project: Contagion2026 / IDMs  
Purpose: GitHub publication refactor without changing scientific computation

---

## 0. Governing principle

The purpose of this work is to move the already frozen scientific workflow into a public, reusable, and reproducible repository structure.

The refactor must not change:

- model definitions
- likelihoods
- parameter transforms
- Gauss-Hermite quadrature rules
- `GH_N`
- optimizer choices
- optimizer tolerances that affect scientific results
- multistart grids
- random seeds
- boundary handling
- clipping conventions
- canonical numerical results

The refactor is limited to:

- moving exact reusable implementations into modules
- simplifying notebook imports
- standardizing paths
- separating public and private data
- improving documentation
- preserving reproducible environments
- adding regression checks
- preparing the repository for public release

If a change appears to require modification of a scientific formula or numerical convention, stop and report it before modifying the code.

---

## 1. Frozen Git baseline

Canonical pre-refactor scientific baseline:

- commit: `ec09459`
- tag: `canonical-pre-refactor-20260906`
- protected branch: `main`
- active refactor branch: `code-publication-refactor`

The baseline contains the complete canonical notebook workflow, Python source files, Julia source/tests, public synthetic data, environment specification, and approved manuscript figures.

Private Moody's data are excluded from Git history.

---

## 2. Current repository structure

Canonical notebooks:

- `01_generate_synthetic_data.ipynb`
- `02_numerical_experiments.ipynb`
- `03A_static_one_period_models.ipynb`
- `03B_environmental_models.ipynb`
- `03C_M3_python.ipynb`
- `03D_M3_julia.ipynb`
- `03E_M3_cross_language_validation.ipynb`
- `03F_subperiod_analysis.ipynb`

Current Python source:

- `models.py`
- `analysis_helpers.py`
- `revision_models.py`

Current module directory:

- `modules/` — currently empty

Julia package:

- `julia/src/Quadrature.jl`
- `julia/src/ProbitNormal.jl`
- `julia/src/Vasicek.jl`
- `julia/src/GaussianMixture.jl`
- `julia/src/IDMLikelihoods.jl`
- scripts under `julia/scripts/`
- tests under `julia/test/`

Data:

- private canonical data: `data/M_1920_2023.csv`
- public synthetic data: `data/M_1920_2023_synthetic.csv`

Outputs:

- manuscript figures: `figure/`
- numerical intermediates and detailed empirical outputs: `pdata/`

Documentation:

- `docs/audits/`
- `docs/implementations/`

---

## 3. Current Python dependency graph

The root modules currently form the following dependency structure:

```text
models.py
   |
   v
analysis_helpers.py
   |
   v
revision_models.py
```

More precisely:

- `analysis_helpers.py` imports functions from `models.py`.
- `revision_models.py` imports `models as m`.
- `revision_models.py` imports `analysis_helpers as ah`.

Observed notebook dependency:

- `03C_M3_python.ipynb` explicitly imports `revision_models.py`.
- `03E_M3_cross_language_validation.ipynb` explicitly imports `revision_models.py`.
- `03A`, `03B`, and `03F` contain substantial inline copies of functions already present in root modules or repeated across notebooks.
- `02` contains numerical utility and model functions duplicated with `models.py`.

The current dependency direction should be preserved during the first-stage refactor.

No circular dependency should be introduced.

Audit records:

- `docs/audits/20260906_python_dependency_and_duplication_audit.txt`
- `docs/audits/20260906_root_module_dependencies.txt`
- `docs/audits/20260906_notebook_module_dependencies.txt`

---

## 4. Julia dependency graph

Julia is already substantially modularized.

Entry point:

```text
IDMLikelihoods.jl
    |
    +-- Quadrature.jl
    +-- ProbitNormal.jl
    +-- Vasicek.jl
    +-- GaussianMixture.jl
```

Important internal dependencies:

- `ProbitNormal.jl` uses `Quadrature`.
- `Vasicek.jl` uses `ProbitNormal.log_comb`.
- `GaussianMixture.jl` uses `Quadrature` and `ProbitNormal.log_comb`.

Scripts and tests load:

```text
julia/src/IDMLikelihoods.jl
```

Canonical Julia environment:

- Julia 1.12.7
- `Project.toml`
- `Manifest.toml`
- 57 / 57 tests PASS
- Phase J2: `GH_N = 100`
- optimizer: Nelder-Mead
- max iterations: 6000
- 180 starts per class
- ALL / SG / IG = 540 starts total

No scientific restructuring of the Julia likelihood layer is planned in this refactor.

Only path handling, documentation, notebook integration, and final reproducibility checks should be changed unless a concrete defect is found.

---

## 5. Duplicate-function audit

The pre-refactor AST audit distinguishes:

- `EXACT-AST`: structurally identical implementations
- `DIFFERENT`: same name but different implementation

This distinction is mandatory during migration.

### 5.1 High-confidence exact duplicates whose existing root implementation is a strong canonical-source candidate

From `models.py`:

- `es_from_pmf`
- `var_from_pmf`
- `moments_from_pmf`
- `m_rho_from_pmf`
- `b_of_s`
- `bisect_root`
- `find_params_on_iso_Torri`
- `pmf_Torri`
- `default_corr_from_asset_corr`
- `invert_asset_corr_from_default_corr`
- `log_pmf_LD_1pt`
- `log_pmf_Torri`
- `pmf_Vas`
- `tail_metrics_from_pmf`

From `analysis_helpers.py`:

- `summarize_nt_default`
- `summarize_nt_default_by_period`
- `fit_one_Torri`
- `neg_ll_torri`
- `fit_one_LD`
- `neg_ll_ld`
- `sigmoid`
- `softplus`
- `inv_softplus`
- `_clean_nh`
- `_clean_nL`
- `gh_normal_weights`
- `scaled_counts`
- `torri_m_q_for_n`
- `ld_m_q_for_n`
- `torri_m_q_fn`
- `ld_m_q_fn`
- `vas_m_q_fn`
- `to_float_dict`

These are the lowest-risk initial candidates for replacing notebook-inline copies with imports.

Function bodies should not be rewritten during that move.

### 5.2 Exact environmental functions duplicated between 03B and 03F

The following are exact duplicates between the environmental and subperiod notebooks:

- `fit_one_pn_LD`
- `fit_one_pn_Torri`
- `fit_pn_null_from_grouped`
- `log_pmf_pn_LD_1pt`
- `log_pmf_pn_Torri_1pt`
- `neg_ll_pn_LD`
- `neg_ll_pn_Torri`
- `neg_ll_pn_null`
- `variance_decomp_pn_ld_fixed_nbar`
- `variance_decomp_pn_torri_fixed_nbar`

Before moving them, compare directly against the current canonical implementation in `revision_models.py`.

If numerically and structurally identical, prefer `revision_models.py` as the canonical source because 03C and 03E already depend on it.

If not identical, do not merge until the scientific reason for the difference is understood.

### 5.3 Same-name DIFFERENT functions requiring manual source inspection

Do not merge automatically:

- `aggregate_m_q_to_rho`
- `build_compare_decomp_table`
- `compute_nbar_for_period`
- `data_m_q_hat`
- `empirical_pmf_from_counts`
- `fit_one_Vas`
- `log_comb`
- `log_pmf_Vas_1pt`
- `logsubexp`
- `logsumexp`
- `logsumexp_arr`
- `logsumexp_vec`
- `m_rho_LD`
- `make_grouped_nh`
- `model_aggregate_m_q`
- `model_pmf_fixed_n`
- `moments_from_params_Torri`
- `neg_ll_vas`
- `normalize_logpmf`
- `phi2_aa`
- `pmf_LD`
- `prepare_plot_data`
- `run_period_analysis`
- `scaled_counts`
- `solve_pq_for_mrho_LD`
- `subset_df_by_period`
- `survival_from_pmf`

For every same-name DIFFERENT case, inspect at least:

1. function signature
2. defaults
3. clipping conventions
4. `GH_N`
5. boundary behavior
6. handling of zero / invalid counts
7. optimizer settings
8. return values
9. fixed-`n` versus observed-`n_t` assumptions
10. downstream tables, figures, or regression audits that use the function

A same-name function may intentionally encode a scientifically different calculation.

### 5.4 Local names that should not be treated as module duplicates

Names such as:

- `F`
- `_fn`
- `objective`
- `residual`
- `unpack`

are local or nested helpers in several contexts.

They must not be promoted to shared global functions merely because the names repeat.

---

## 6. Proposed Python module structure

The first publication refactor should be conservative.

Proposed structure:

```text
modules/
├── __init__.py
├── models.py
├── analysis_helpers.py
├── revision_models.py
└── data_io.py
```

The first migration should preserve the existing three-layer Python architecture rather than redesign it.

### 6.1 `modules/models.py`

Initial source:

- exact copy of the current `models.py`

Responsibilities:

- one-period Torri functions
- Davis-Lo functions
- Vasicek functions
- PMF / survival functions
- moment utilities
- VaR / ES utilities
- iso-(m,rho) calculations
- low-level numerical likelihood utilities

No scientific body should be rewritten during the first migration.

### 6.2 `modules/analysis_helpers.py`

Initial source:

- exact copy of the current `analysis_helpers.py`

Only import paths should initially change, for example:

```python
from .models import ...
```

Responsibilities:

- fitting wrappers
- grouped-data helpers
- descriptive statistics
- variance-decomposition support
- plotting preparation
- period helpers
- model-comparison table helpers

### 6.3 `modules/revision_models.py`

Initial source:

- exact copy of the current `revision_models.py`

Only import paths should initially change:

```python
from . import models as m
from . import analysis_helpers as ah
```

Responsibilities:

- canonical revision likelihoods
- M0/M1/M2/M3 functions
- canonical M3 physical likelihood
- parameter transforms
- M3 multistart support

The filename `revision_models.py` should be retained during this publication refactor to minimize dependency and numerical risk.

Renaming can be considered only after the publication refactor is complete and is not required.

### 6.4 `modules/data_io.py`

New non-scientific utility module.

Responsibilities:

- project-root resolution
- data directories
- figure directory
- `pdata` directory
- explicit real/synthetic selection
- schema validation
- data-mode reporting

The module must not transform scientific values beyond the existing loading conventions.

Suggested constants:

```text
PROJECT_ROOT
DATA_DIR
FIGURE_DIR
PDATA_DIR
```

No absolute user path is allowed.

---

## 7. Root-module migration policy

Do not immediately delete:

- `models.py`
- `analysis_helpers.py`
- `revision_models.py`

Recommended migration:

1. copy the three root modules into `modules/`
2. adjust only their internal import paths
3. compare imported functions with the root versions
4. run import smoke tests
5. migrate notebooks one at a time
6. run notebook-level numerical regression after each migration
7. remove root modules only after all notebooks pass

If necessary, temporary compatibility wrappers may be used, but they should not survive into the final public repository unless justified.

---

## 8. Notebook migration plan

### Stage target 01 — `01_generate_synthetic_data.ipynb`

Scope:

- path handling only if needed
- optional use of `modules.data_io`
- preserve generation algorithm
- preserve random seed
- preserve public output filename:

```text
data/M_1920_2023_synthetic.csv
```

No statistical change.

### Stage target 02 — `02_numerical_experiments.ipynb`

This is the first scientific notebook migration.

Replace only confirmed EXACT-AST inline functions with imports from:

- `modules.models`
- `modules.analysis_helpers`

Do not consolidate DIFFERENT implementations.

Special caution:

- `F`
- `m_rho_LD`
- `pmf_LD`
- `solve_pq_for_mrho_LD`
- `objective`
- `survival_from_pmf`

Run the notebook top-to-bottom and compare all numerical outputs and generated figures against the baseline.

### Stage target 03A — `03A_static_one_period_models.ipynb`

Replace exact static-model and descriptive helper duplicates first.

Manual inspection required before replacing:

- `log_comb`
- `logsumexp`
- `logsubexp`
- `normalize_logpmf`
- `fit_one_Vas`
- `log_pmf_Vas_1pt`
- `neg_ll_vas`
- aggregate moment/correlation helpers

Table I, Table II, and static benchmark figures must remain unchanged.

### Stage target 03B — `03B_environmental_models.ipynb`

Replace only exact duplicates whose canonical source is unambiguous.

Then compare environmental M0/M1/M2-related functions against `revision_models.py`.

No environmental likelihood is to be changed merely for code cleanup.

Verify all full-sample regression targets and boundary audits.

### Stage target 03C — `03C_M3_python.ipynb`

This notebook already treats `revision_models.py` as canonical.

Migration should be minimal:

```python
import modules.revision_models as rm
```

Preserve:

- physical M3 likelihood
- GH_N
- start generation
- multistart optimizer
- boundary logic
- reported Python local best

### Stage target 03D — `03D_M3_julia.ipynb`

No scientific Julia refactor.

Only verify:

- relative path handling
- explicit real-data requirement
- public/private data messaging
- invocation of the current Julia project/scripts

The 540-start Phase J2 configuration must remain unchanged.

### Stage target 03E — `03E_M3_cross_language_validation.ipynb`

Migration should be minimal:

```python
import modules.revision_models as rm
```

Preserve:

- canonical Python physical likelihood
- Python/Julia cross-evaluation
- selected candidate logic
- final source selection

Required terminology:

```text
cross-language validated best-found candidate
```

Do not use:

```text
global optimum
```

### Stage target 03F — `03F_subperiod_analysis.ipynb`

Migrate only after 03A and 03B are stable.

It contains substantial exact duplication with 03B and `analysis_helpers.py`.

Subperiod-specific DIFFERENT helpers require manual source comparison before consolidation.

Preserve:

- subperiod NLLs
- boundary classifications
- variance decomposition
- PMF checks
- canonical regression audits

---

## 9. Data-path policy

All public code must use repository-relative paths.

Private canonical data:

```text
data/M_1920_2023.csv
```

Public synthetic data:

```text
data/M_1920_2023_synthetic.csv
```

The private file remains Git-ignored.

No code, documentation, or notebook output may contain:

```text
[LOCAL_PATH]
Dropbox/
```

### 9.1 Real/synthetic behavior

The public repository contains only synthetic data.

Synthetic data do not reproduce manuscript numerical estimates.

For workflows where synthetic execution is scientifically meaningful, explicit selection may be used.

Canonical manuscript audits must never silently substitute synthetic data for real data.

In particular:

- 03D canonical Julia Phase J2 requires real data
- 03E canonical cross-language manuscript validation requires real data
- manuscript-specific regression audits must report their data mode explicitly

If synthetic demonstration mode is allowed, print an explicit message such as:

```text
DATA MODE: synthetic demonstration
Paper numerical estimates will not be reproduced.
```

Silent fallback is prohibited.

---

## 10. Output-path policy

Standard directories:

```text
figure/
docs/
data/
pdata/
modules/
julia/
```

Do not create `figures/`.

Canonical figures remain in `figure/`.

Examples:

- `fig_survival_hierarchical_ALL_IG.pdf`
- `fig_pmf_hier_ALLCLASS.pdf`
- `fig_variance_decomp_hier_compare.pdf`
- `fig_variance_decomp_two_periods.pdf`

Detailed generated numeric outputs continue to use `pdata/`.

`pdata/` remains excluded from Git unless individual files are separately reviewed and approved.

---

## 11. Public/private data separation

### Public

- Python source
- notebooks
- Julia source/scripts/tests
- synthetic dataset
- approved manuscript figures
- reviewed compact regression reference values
- environment files
- documentation

### Private

- `data/M_1920_2023.csv`
- raw proprietary Moody's observations
- unreviewed detailed empirical outputs under `pdata/`

### Approved manuscript-derived public artifacts

The current manuscript figures under `figure/` have been approved for release.

`julia/test/reference_python_results.csv` contains only compact fixed-parameter regression values and is approved as a public numerical test fixture.

---

## 12. `.gitignore` policy

Required exclusions include:

```text
.DS_Store
__pycache__/
*.py[cod]
.pytest_cache/
.mypy_cache/
.ruff_cache/
.ipynb_checkpoints/
.venv/
venv/
.env
.env.*
data/M_1920_2023.csv
pdata/
*.tmp
*.bak
*~
```

Before any public push, verify both the Git index and Git history.

---

## 13. README plan

Create a root `README.md` containing:

1. paper/project title
2. project purpose
3. repository structure
4. notebook execution order
5. Python environment
6. Julia 1.12.7 environment
7. synthetic-data policy
8. proprietary-data exclusion
9. exact manuscript reproduction requirements
10. output directories
11. M3 Python/Julia cross-language validation
12. meaning of "cross-language validated best-found candidate"
13. licensing placeholder
14. citation placeholder

Required statement:

> The public repository contains synthetic data only. The synthetic dataset is provided to demonstrate the code and workflow and does not reproduce the numerical estimates reported in the paper. Exact reproduction of manuscript estimates requires the original same-schema proprietary dataset, which is not included in this repository.

Do not describe how to obtain proprietary data.

---

## 14. Python environment reproducibility

Current environment:

- `environment.yml`
- local absolute `prefix:` already removed

Policy:

- do not upgrade packages during this refactor
- do not replace package versions because newer versions exist
- preserve the working canonical environment
- add missing explicit dependencies only if required for clean recreation

After refactor:

1. create a fresh environment from `environment.yml`
2. run import smoke tests
3. execute notebooks in documented order
4. run numerical regression checks

A `requirements.txt` should only be added if it materially improves portability without altering the working dependency set.

---

## 15. Canonical numerical regression targets

All targets below must survive refactoring.

### 15.1 03B full sample

Table III:

```text
ALL
M0 = 431.363710578
M1 ≈ 413.825
M2 = 431.363710578

SG
M0 = 417.708534945
M1 ≈ 401.342
M2 = 409.266628263

IG
M0 = 182.381686046
M1 ≈ 180.781
M2 = 182.381686046
```

Boundary requirements:

```text
ALL M2 = M0 at v = 0
IG  M2 = M0 at v = 0
SG  M2 is an interior solution
```

At the ALL and IG M2 boundary:

- `u` is scientifically unidentified
- PMF equality with M0 must be retained
- survival equality with M0 must be retained
- derived risk-measure equality must be retained where tested

Table VII M1:

```text
ALL .026 / .090 / .575
SG  .034 / .103 / .648
IG  .137 / .101 / 1.175
```

Table VII M2:

```text
ALL .025 / .000 / .998
SG  .041 / .019 / 1.666
IG  .145 / .000 / 1.826
```

### 15.2 M3

Final best-found NLL:

```text
ALL 430.271292...
SG  416.177207...
IG  180.112948...
```

Final selected source:

```text
ALL: Julia
SG:  Julia
IG:  Python
```

Python SG local best:

```text
416.263555...
```

This is intentionally worse than the Julia candidate and must not be changed to match Julia.

### 15.3 Subperiods

1950-1979:

```text
ALL M0/M2 58.2442
ALL M1    58.1637

SG M0     58.8513
SG M1     58.8388
SG M2     55.9583

IG M0     12.9056
IG M1     11.8125
IG M2     11.7179
```

1980-2023:

```text
ALL M0/M2 214.613
ALL M1    212.726

SG M0/M1/M2 208.747

IG M0/M2 84.7066
IG M1    83.4159
```

### 15.4 Julia

Tests:

```text
57 / 57 PASS
```

Canonical Phase J2:

```text
ALL 430.2712920847
SG  416.1772074557
IG  180.1129479483
```

---

## 16. Regression methodology

Every migration stage must follow:

1. record current Git commit
2. modify only the targeted module/notebook
3. inspect `git diff`
4. run import smoke tests
5. execute the targeted notebook
6. compare canonical tables and diagnostics
7. compare boundary classifications
8. verify generated output filenames
9. commit only after regression passes

For EXACT-AST duplicates, source movement/import substitution is preferred over rewriting.

For DIFFERENT duplicates, perform source diff before deciding whether to consolidate.

Do not apply broad formatters or cleanup tools to scientific source until numerical regression is complete.

---

## 17. Refactor stages

### Stage 0 — completed

- Git repository initialized
- private data ignored
- privacy/path audit performed
- absolute local paths removed from tracked content
- canonical pre-refactor baseline committed
- canonical baseline tagged
- refactor branch created
- dependency/duplication audits generated

### Stage 1 — module scaffolding

Create:

- `modules/__init__.py`
- `modules/models.py`
- `modules/analysis_helpers.py`
- `modules/revision_models.py`
- `modules/data_io.py`

Initially preserve existing scientific function bodies exactly.

Run import smoke tests.

Do not modify notebooks in the same commit.

### Stage 2 — Notebook 02 migration

Replace only confirmed exact duplicates.

Run 02 top-to-bottom.

Compare all numerical outputs and figures.

### Stage 3 — Notebook 03A migration

Migrate exact static-model helpers.

Manually inspect DIFFERENT Vasicek and low-level numerical functions.

Run full regression.

### Stage 4 — Notebook 03B migration

Migrate exact environmental functions only after comparison with canonical `revision_models.py`.

Verify Table III, Table IV, Table VII, Table VIII, boundary equality, PMF/survival diagnostics.

### Stage 5 — 03C / 03E migration

Change canonical import to `modules.revision_models`.

Verify M3 likelihood, cross-evaluation, and candidate selection.

Do not alter optimization logic.

### Stage 6 — 03F migration

Migrate shared helpers only after 03B is stable.

Run all subperiod, boundary, and regression audits.

### Stage 7 — Julia integration

Run:

- Julia unit tests
- fixed-parameter validation
- cross-language validation

Do not change Julia likelihood logic.

The full 540-start Phase J2 reproduction may be deferred until the final complete run unless a Julia scientific source file changes.

### Stage 8 — Notebook 01 and complete fresh run

Run the documented public workflow from a clean environment.

Confirm:

- no private path leakage
- no silent synthetic fallback in manuscript audits
- outputs use documented directories
- all regression targets pass

---

## 18. Git strategy

Use small stage-specific commits.

Suggested commit sequence:

```text
Add pre-refactor dependency and duplication audits
Add code publication refactor implementation plan
Add Python module scaffold
Migrate numerical experiments notebook to shared modules
Migrate static benchmark notebook to shared modules
Migrate environmental notebook to shared modules
Migrate M3 Python and cross-language imports
Migrate subperiod notebook to shared modules
Standardize public data-path handling
Add repository README and reproduction guide
```

Do not mix scientific changes with formatting changes.

Whitespace cleanup should be separate and preferably occur only after scientific regression is complete.

---

## 19. GitHub pre-publication audit

Before any push, search the working tree, index, and relevant Git history for:

- `M_1920_2023.csv`
- proprietary Moody's data
- raw empirical records
- `[LOCAL_PATH]`
- `Dropbox`
- credentials
- passwords
- API keys
- tokens
- `.DS_Store`
- `__pycache__`
- temporary files
- private notebook output
- stale private paths
- unreviewed `pdata` files

Verify explicitly that:

```text
data/M_1920_2023.csv
```

has never entered tracked Git history.

No GitHub push is permitted until this audit passes and explicit approval is given.

---

## 20. Expected final public repository structure

```text
.
├── README.md
├── environment.yml
├── 01_generate_synthetic_data.ipynb
├── 02_numerical_experiments.ipynb
├── 03A_static_one_period_models.ipynb
├── 03B_environmental_models.ipynb
├── 03C_M3_python.ipynb
├── 03D_M3_julia.ipynb
├── 03E_M3_cross_language_validation.ipynb
├── 03F_subperiod_analysis.ipynb
├── data/
│   └── M_1920_2023_synthetic.csv
├── figure/
├── modules/
│   ├── __init__.py
│   ├── models.py
│   ├── analysis_helpers.py
│   ├── revision_models.py
│   └── data_io.py
├── julia/
│   ├── Project.toml
│   ├── Manifest.toml
│   ├── README.md
│   ├── src/
│   ├── scripts/
│   └── test/
└── docs/
    ├── audits/
    └── implementations/
```

`pdata/` may remain present locally but remains excluded from Git.

---

## 21. Stop conditions

Stop the refactor and report before proceeding if any of the following occurs:

- a canonical NLL changes beyond the established numerical tolerance
- a boundary solution changes classification
- M0/M2 equality fails for ALL or IG
- Python/Julia fixed-parameter likelihoods diverge
- SG Python unexpectedly changes relative to the canonical Julia J2 candidate
- `GH_N` changes
- optimizer settings change
- clipping rules change
- parameter transforms change
- a notebook requires a scientific-formula change to complete migration
- proprietary data enter Git
- synthetic data are silently used for a manuscript audit

---

## 22. Governing completion criterion

The refactor is complete only when:

> The scientific computation is unchanged relative to `canonical-pre-refactor-20260906`, reusable functions are imported from a publication-ready module structure, private data remain outside Git, and the complete documented workflow passes its numerical regression checks.
