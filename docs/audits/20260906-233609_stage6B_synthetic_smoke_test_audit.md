# Stage 6B Audit: Public Synthetic-State Smoke Test

## Stage 6B v1 Finding
- **Execution**: The synthetic workflow was executed top-to-bottom (01, 02, 03A, 03B, 03C, 03F).
- **Result**: 01 through 03C PASSED. `03F_subperiod_analysis.ipynb` FAILED during boundary assertion (`assert bool(boundary_audit["pass"].all())`).
- **Cause**: Although 03F has a synthetic fallback loader, the downstream boundary audit assertion mathematically requires the boundary structure specifically obtained from the proprietary real data. The synthetic data does not replicate this boundary structure.
- **Action**: 03F was removed from the README's supported public synthetic demonstration workflow. No scientific source code or assertions were modified. The Stage 6A real-data PASS is unaffected.

## Stage 6B v2 Validation
- **Public Tree Creation**: A clean, Git-tracked-only archive was extracted to `/tmp/idms_public_smoke_v2`.
- **Privacy Check**: Verified that `data/M_1920_2023.csv` and `pdata/` do NOT exist in the public archive. `data/M_1920_2023_synthetic.csv` is correctly present.
- **Workflow Execution**: The updated validated synthetic workflow (01, 02, 03A, 03B, 03C) was executed sequentially using the `idms-jpsj` kernel.
- **Result**: All notebooks executed successfully (exit code 0) without referring to proprietary data. The outputs explicitly displayed "FULLY SYNTHETIC" warnings as expected.
- **Source Integrity**: Code and markdown cells remained unchanged.

## Public-Policy Verification
- `03D_M3_julia.ipynb` explicitly requires real data and does not claim synthetic manuscript validation.
- `03E_M3_cross_language_validation.ipynb` relies on outputs generated from real data and does not claim synthetic validation.
- `03F_subperiod_analysis.ipynb` is no longer claimed as a validated synthetic workflow in the README.

**Result**: Stage 6B PASS
