# Stage 5A: Publication Readiness Audit (Final)

## 1. Git/Environment Preflight
- **Branch**: `code-publication-refactor`
- **Status**: Clean
- **Canonical Python**: `idms-jpsj` env.
- **Environment Invocation Policy**: To ensure robust, path-independent reproducibility, all future environment invocations must use `conda run -n idms-jpsj` instead of absolute executable paths. Example:
  `conda run -n idms-jpsj jupyter nbconvert --ExecutePreprocessor.kernel_name=idms-jpsj ...`
- **Environment Restoration**: The `pyyaml` package was inadvertently installed during early audit steps. It has been strictly uninstalled (`pip uninstall -y pyyaml`), restoring the environment to its exact canonical lock state. Verification (`import numpy, scipy, pandas, matplotlib`) passed.

## 2. Local Absolute Paths Audit
Found 15 hits of `/Users/mori` in tracked files:
- **[B. notebook output metadata]**: `02_numerical_experiments.ipynb`, `03A_static_one_period_models.ipynb` (e.g. `figure output: ...`, `Saved: .../pdata/...`)
- **[C. audit/report provenance]**: Past audit markdown logs containing local path references.
- **[D. harmless textual example]**: Textual examples in implementation plans.
- **Action**: These are **REQUIRED BEFORE PUSH** (must remove before publication) to sanitize the repository, though they pose no execution portability or security risk.

## 3. Git History Privacy Audit
- **Proprietary Data (`data/M_1920_2023.csv`)**: ZERO hits in Git history. Secure.
- **Local Paths (`/Users/mori`, `Dropbox`)**: 15 hits exist in both the current tracked tree AND in Git history commits.
- **`pdata/`**: ZERO hits in Git history.
- **`.DS_Store` / `__pycache__`**: ZERO hits in Git history.

## 4. Notebook Data Mode Table
| Notebook | Real Data Req? | Synth Supported? | Fallback Exists? | Fallback Type | Manuscript Reproducible w/ Synth? | pdata In | pdata Out |
|---|---|---|---|---|---|---|---|
| 01 | No | Yes | No (generates) | N/A | N/A | No | No |
| 02 | No | N/A | No | N/A | Yes (independent) | No | Yes |
| 03A | Yes | Yes | Yes | Explicit Warning | No | Yes | Yes |
| 03B | Yes | Yes | Yes | Explicit Warning | No | Yes | Yes |
| 03C | Yes | Yes | Yes | Explicit Warning | No | Yes | Yes |
| 03D | Yes | No | No | N/A | No | No | No |
| 03E | Yes | No | No | N/A | No | Yes | Yes |
| 03F | Yes | Yes | Yes | Explicit Warning | No | Yes | Yes |

*Note: 03F is completely independent of 03B and does not read any `pdata`. 03D (Julia M3) and 03E (Cross-language Validation) strictly require real data and do not support synthetic fallback.*

## 5. Execution Workflows
**A. Exact Manuscript Reproduction Workflow**
*(Requires proprietary `data/M_1920_2023.csv`)*
- `02_numerical_experiments.ipynb` (independent)
- `03A_static_one_period_models.ipynb` (independent)
- `03B_environmental_models.ipynb` (independent, produces downstream `pdata`)
- `03C_M3_python.ipynb` (depends on `03B` outputs)
- `03D_M3_julia.ipynb` (independent Julia M3 manuscript run using real data)
- `03E_M3_cross_language_validation.ipynb` (depends on outputs/results needed from `03B/03C/03D`)
- `03F_subperiod_analysis.ipynb` (independent subperiod analysis using real data directly)
*(Note: `01_generate_synthetic_data.ipynb` is NOT required for exact manuscript reproduction).*

**B. Public Synthetic Demonstration Workflow**
*(Requires no real data)*
- `01_generate_synthetic_data.ipynb` (generates synthetic data)
- `02_numerical_experiments.ipynb`
- `03A_static_one_period_models.ipynb`
- `03B_environmental_models.ipynb`
- `03C_M3_python.ipynb`
- `03F_subperiod_analysis.ipynb`
*(Note: `03D` and `03E` are EXCLUDED from the synthetic manuscript-validation workflow).*

## 6. Python Environment Audit
| Package | Status in `environment.yml` | Usage |
|---|---|---|
| `numpy`, `scipy`, `pandas`, `matplotlib`, `ipykernel`, `python=3.11` | Explicitly listed | Direct Import / Core runtime |
| `nbconvert` / `jupyter` | **NOT EXPLICITLY LISTED** | Available only transitively |

*Note: `nbconvert` is missing as a direct dependency. While functional in the current environment, fresh environment recreation is not guaranteed. **Stage 5B should evaluate adding `nbconvert` as a direct dependency**.*

## 7. Root vs Modules Assessment
- **State**: `02` imports `modules.models`. `03A` imports `modules.analysis_helpers`. `03C`/`03E` import **root** `revision_models`.
- **Conclusion**: Root files (`revision_models.py`, `models.py`, `analysis_helpers.py`) MUST be maintained. Deleting them breaks `03C` and `03E` and destroys exact reproducibility.

## 8. Blocker Classification

**CRITICAL (0 issues):**
- Proprietary data leakage: None.
- Credentials/secrets: None.

**REQUIRED BEFORE PUSH (4 issues):**
- Tracked local absolute paths cleanup (`[LOCAL_REPO_PATH]/...`).
- Inconsistent README/data-mode description (README missing).
- Missing final full reproduction run (using real data top-to-bottom).
- Missing final public-state smoke test (synthetic workflow run).

**OPTIONAL / NON-BLOCKING:**
- `LICENSE` / `CITATION.cff` metadata (highly recommended for release).
- Cosmetic duplication (root vs modules, duplicated PNG vs PDF).
- Further scientific-code modularization.
