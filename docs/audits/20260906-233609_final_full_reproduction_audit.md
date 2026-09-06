# Final Full Real-Data Reproduction Audit

## 1. Environment & Preflight
- **Branch**: `code-publication-refactor`
- **Starting Commit**: 8834f06
- **Python Environment**: `idms-jpsj` (Python 3.11.16)
- **External Runner**: `jupyter-nbconvert` (7.17.1) at `/Library/Frameworks/Python.framework/Versions/3.11/bin/jupyter-nbconvert`
- **Kernel Confirmation**: Notebooks 01-03C, 03E, 03F successfully utilized `--ExecutePreprocessor.kernel_name=idms-jpsj`. Notebook 03D successfully utilized `--ExecutePreprocessor.kernel_name=julia-1.12`.
- **Julia Version**: 1.12.7
- **Real-Data Mode**: Confirmed active (no synthetic fallback warnings during 03A, 03B, 03C, 03E, 03F).

## 2. Execution Status per Component
- **01_generate_synthetic_data**: PASS. Synthetic CSV hash unchanged.
- **02_numerical_experiments**: PASS. 
- **03A_static_one_period**: PASS. 
- **03B_environmental_models**: PASS. M0/M1/M2 boundary checks conform to exact targets (e.g., M0=431.363 ALL, boundary equivalent to M2).
- **03C_M3_python**: PASS. Correct SG local best ~416.263 preserved.
- **Julia Test Suite**: PASS. 57/57 tests.
- **03D_M3_julia**: PASS. 540 starts total (180 per class). ALL ~430.271, SG ~416.177, IG ~180.113. 
- **03E_M3_cross_language_validation**: PASS. Selected Julia for ALL/SG, Python for IG.
- **03F_subperiod_analysis**: PASS. Correct subperiod target regressions.

## 3. Source Integrity & Numerical Regression
- **Code & Markdown Source**: 100% UNCHANGED for all notebooks (pre-run hashes exactly match post-run hashes).
- **Numerical Regression**: All sanity targets strictly met (abs diff <= 1e-6).
- **M3 Terminology**: M3 identified properly as the cross-language validated best-found candidate.

## 4. Tracked File Handling
- Tracked output changes (notebook metadata, re-generated PDF figures) were generated cleanly but subsequently restored to HEAD to preserve execution-only metadata stability.
- No proprietary data or secrets tracked.

## Overall Result
**PASS**
