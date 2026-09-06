# IDMLikelihoods.jl - Phase J1 & J2 Cross-Language Validation

This package provides the Julia high-performance numerical engine for credit default contagion models, supporting cross-language validation against Python canonical implementations. The canonical Julia version for this project is **1.12.7**.

## Phase J1 Scope
- Single Probit-Normal NULL model (**M0**, $k=2$).
- Single-factor **Vasicek** model.
- Bi-directional parameter transformation:
  $$\mu = \frac{\Phi^{-1}(p_0)}{\sqrt{1 - \rho_A}}, \quad \sigma = \sqrt{\frac{\rho_A}{1 - \rho_A}}$$
  $$p_0 = \Phi\left(\frac{\mu}{\sqrt{1 + \sigma^2}}\right), \quad \rho_A = \frac{\sigma^2}{1 + \sigma^2}$$

## Phase J2 Scope / Gaussian Mixture M3
- Supports advanced Gaussian mixture modeling for contagion components (**M3**).
- Includes scripts for fitting mixtures (`fit_mixture.jl`) and testing the mixture formulations (`test_gaussian_mixture.jl`).
- Plays a critical cross-language validation role to independently optimize and cross-evaluate parameters with Python, verifying the **cross-language validated best-found candidate** rather than guaranteeing a theoretical global optimum.

## Running Cross-Language Validation Scripts

1. **Fixed-Parameter Validation**:
   ```bash
   julia --project=julia julia/scripts/validate_fixed_parameters.jl
   ```

2. **M0 MLE Optimization & Benchmark Export**:
   ```bash
   julia --project=julia julia/scripts/fit_m0_mle.jl
   ```

3. **Automated Unit & Integration Test Suite**:
   ```bash
   julia --project=julia -e 'using Pkg; Pkg.test()'
   ```
   *(Note: The canonical test target verifies that 57 / 57 tests PASS).*
