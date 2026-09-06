import math
import builtins
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp as scipy_logsumexp
from numpy.polynomial.hermite import hermgauss

def logsumexp_two_scalars(a, b):
    m = max(a, b)
    if math.isinf(m) and m < 0:
        return -math.inf
    return m + math.log(math.exp(a - m) + math.exp(b - m))

def scalar_or_array_logsumexp(a, b=None, **kwargs):
    if b is not None and isinstance(b, (int, float, np.number)):
        return logsumexp_two_scalars(float(a), float(b))
    return scipy_logsumexp(a, b=b, **kwargs)

builtins.logsumexp = scalar_or_array_logsumexp

def logsumexp_vec(logw):
    m_val = np.max(logw)
    if np.isneginf(m_val):
        return -np.inf
    return float(m_val + np.log(np.sum(np.exp(logw - m_val))))

builtins.logsumexp_vec = logsumexp_vec

from . import models as m
from . import analysis_helpers as ah

def vasicek_to_probit_normal(p0: float, rho_A: float):
    """
    Convert Vasicek parameters (p0, rho_A) to Probit-Normal (mu, sigma) parameters.
    mu = norm.ppf(p0) / sqrt(1 - rho_A)
    sigma = sqrt(rho_A / (1 - rho_A))
    """
    if not (0.0 < p0 < 1.0) or not (0.0 <= rho_A < 1.0):
        raise ValueError("Require 0 < p0 < 1 and 0 <= rho_A < 1")
    
    val_inv = math.sqrt(1.0 - rho_A)
    mu = norm.ppf(p0) / val_inv
    sigma = math.sqrt(rho_A / (1.0 - rho_A))
    return float(mu), float(sigma)

def probit_normal_to_vasicek(mu: float, sigma: float):
    """
    Convert Probit-Normal parameters (mu, sigma) to Vasicek (p0, rho_A) parameters.
    p0 = norm.cdf(mu / sqrt(1 + sigma^2))
    rho_A = sigma^2 / (1 + sigma^2)
    """
    if sigma < 0:
        raise ValueError("sigma must be >= 0")
    
    var_y = sigma**2
    p0 = norm.cdf(mu / math.sqrt(1.0 + var_y))
    rho_A = var_y / (1.0 + var_y)
    return float(p0), float(rho_A)

def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -float("inf")
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))

def nll_probit_normal(params, n_arr, h_arr, c_arr, GH_N=100):
    mu, log_sigma = params
    sigma = math.exp(log_sigma)

    x, w = hermgauss(int(GH_N))
    y_nodes = mu + math.sqrt(2.0) * sigma * x
    w_nodes = w / math.sqrt(math.pi)
    w_nodes = w_nodes / np.sum(w_nodes)

    p_nodes = norm.cdf(y_nodes)
    p_nodes = np.clip(p_nodes, 1e-15, 1.0 - 1e-15)

    log_p = np.log(p_nodes)
    log_1mp = np.log1p(-p_nodes)
    log_w = np.log(w_nodes)

    tot_nll = 0.0
    for n, h, cnt in zip(n_arr, h_arr, c_arr):
        log_c = log_comb(int(n), int(h))
        log_terms = log_c + h * log_p + (n - h) * log_1mp
        log_prob = scalar_or_array_logsumexp(log_w + log_terms)
        tot_nll -= cnt * log_prob

    return tot_nll

def fit_probit_normal_null(n_arr, h_arr, c_arr=None, GH_N=100):
    """
    Fit Probit-Normal NULL model (M0) via Maximum Likelihood.
    Returns dict with mu, sigma, p0_equiv, rhoA_equiv, nll, success.
    """
    if c_arr is None:
        tmp = pd.DataFrame({"n": n_arr, "h": h_arr})
        g = tmp.value_counts().reset_index(name="cnt")
        n_arr = g["n"].to_numpy(dtype=int)
        h_arr = g["h"].to_numpy(dtype=int)
        c_arr = g["cnt"].to_numpy(dtype=int)

    p_hat = np.sum(h_arr * c_arr) / np.sum(n_arr * c_arr)
    mu_init = float(norm.ppf(p_hat))
    sigma_init = 0.3

    res = minimize(
        nll_probit_normal,
        x0=[mu_init, math.log(sigma_init)],
        args=(n_arr, h_arr, c_arr, GH_N),
        method="L-BFGS-B"
    )

    mu_fit = float(res.x[0])
    sigma_fit = float(math.exp(res.x[1]))
    p0_eq, rhoA_eq = probit_normal_to_vasicek(mu_fit, sigma_fit)

    return {
        "mu": mu_fit,
        "sigma": sigma_fit,
        "p0_equiv": p0_eq,
        "rho_A_equiv": rhoA_eq,
        "nll": float(res.fun),
        "success": bool(res.success),
        "message": str(res.message)
    }

def compute_aic_bic(nll: float, k: int, T: int):
    """
    Compute AIC and BIC from negative log-likelihood (nll), parameter count (k), and actual annual observation count (T = sum(c_arr)).
    """
    logL = -float(nll)
    aic = 2.0 * k - 2.0 * logL
    bic = float(k) * math.log(float(T)) - 2.0 * logL
    return {
        "logL": logL,
        "nll": float(nll),
        "k": int(k),
        "T": int(T),
        "AIC": float(aic),
        "BIC": float(bic)
    }

def test_direct_boundary_nestedness(df_row_series_n, df_row_series_h, GH_N=100):
    """
    Directly evaluate Torri (at v=0) and Davis-Lo (at q=0) using M0 (mu, sigma)
    to confirm exact likelihood equivalence before running multi-start optimization.
    Contagion-free boundary for Torri is v = 0 (or u = 1).
    Contagion-free boundary for Davis-Lo is q = 0.
    """
    n_arr, h_arr, c_arr = ah.make_grouped_nh(df_row_series_n, df_row_series_h)
    fit_m0 = fit_probit_normal_null(n_arr, h_arr, c_arr, GH_N=GH_N)
    mu0, sigma0 = fit_m0["mu"], fit_m0["sigma"]
    nll_m0 = fit_m0["nll"]

    gh_x, gh_w = hermgauss(GH_N)

    # 1. Torri at v=0, u=0.5
    nll_torri_v0 = ah.neg_ll_hier_Torri(
        [mu0, ah.inv_softplus(sigma0), math.log(0.5 / 0.5), math.log(1e-12 / (1.0 - 1e-12))],
        n_arr, h_arr, c_arr, gh_x, gh_w
    )

    # 2. Davis-Lo at q=0
    nll_ld_q0 = ah.neg_ll_hier_LD(
        [mu0, ah.inv_softplus(sigma0), math.log(1e-12 / (1.0 - 1e-12))],
        n_arr, h_arr, c_arr, gh_x, gh_w
    )

    delta_torri = abs(nll_torri_v0 - nll_m0)
    delta_ld = abs(nll_ld_q0 - nll_m0)

    return {
        "nll_m0": nll_m0,
        "nll_torri_v0": nll_torri_v0,
        "nll_ld_q0": nll_ld_q0,
        "delta_torri": delta_torri,
        "delta_ld": delta_ld,
        "torri_nested_ok": bool(delta_torri < 1e-4),
        "ld_nested_ok": bool(delta_ld < 1e-4),
    }

def fit_one_hier_Torri_multistart(series_n, series_h, GH_N=100, maxiter=6000):
    """
    Multi-start optimization for Hierarchical Torri (M2) to ensure finding best-found candidate
    without getting trapped in local minima above M0.
    """
    n_arr, h_arr, c_arr = ah.make_grouped_nh(series_n, series_h)
    gh_x, gh_w = hermgauss(GH_N)

    fit_m0 = fit_probit_normal_null(n_arr, h_arr, c_arr, GH_N=GH_N)
    mu0, sigma0 = fit_m0["mu"], fit_m0["sigma"]

    x0_candidates = []

    # 1. M0 boundary candidates (v -> 0)
    for u0_val in [0.5, 0.9]:
        for v0_val in [1e-6, 1e-3]:
            x0_candidates.append([
                mu0,
                ah.inv_softplus(sigma0),
                math.log(u0_val / (1.0 - u0_val)),
                math.log(v0_val / (1.0 - v0_val))
            ])

    # 2. Paper_V1 crude default
    m_hat = float(series_h.sum() / series_n.sum())
    mu_crude = float(norm.ppf(np.clip(m_hat, 1e-6, 1.0 - 1e-6)))
    x0_candidates.append([
        mu_crude,
        ah.inv_softplus(0.5),
        math.log(0.9 / 0.1),
        math.log(0.01 / 0.99)
    ])

    # 3. Streamlined internal grid
    for u_c in [0.2, 0.7, 0.95]:
        for v_c in [0.01, 0.15]:
            x0_candidates.append([
                mu0,
                ah.inv_softplus(sigma0),
                math.log(u_c / (1.0 - u_c)),
                math.log(v_c / (1.0 - v_c))
            ])

    best_res = None
    best_nll = float("inf")

    for x0 in x0_candidates:
        try:
            res = minimize(
                ah.neg_ll_hier_Torri,
                x0=np.array(x0, dtype=float),
                args=(n_arr, h_arr, c_arr, gh_x, gh_w),
                method="Nelder-Mead",
                options={"maxiter": maxiter}
            )
            if res.fun < best_nll:
                best_nll = float(res.fun)
                best_res = res
        except Exception:
            continue

    mu_hat = float(best_res.x[0])
    sigma_hat = float(ah.softplus(best_res.x[1]))
    u_hat = float(ah.sigmoid(best_res.x[2]))
    v_hat = float(ah.sigmoid(best_res.x[3]))

    return {
        "mu": mu_hat,
        "sigma": sigma_hat,
        "u": u_hat,
        "v": v_hat,
        "nll": best_nll,
        "success": bool(best_res.success),
        "nit": int(best_res.nit),
        "message": str(best_res.message),
        "GH_N": int(GH_N),
    }

def fit_one_hier_LD_multistart(series_n, series_h, GH_N=100, maxiter=6000):
    """
    Multi-start optimization for Hierarchical Davis-Lo (M1) including M0 boundary candidate (q -> 0).
    """
    n_arr, h_arr, c_arr = ah.make_grouped_nh(series_n, series_h)
    gh_x, gh_w = hermgauss(GH_N)

    fit_m0 = fit_probit_normal_null(n_arr, h_arr, c_arr, GH_N=GH_N)
    mu0, sigma0 = fit_m0["mu"], fit_m0["sigma"]

    x0_candidates = []

    # 1. M0 boundary candidates (q -> 0)
    for q0_val in [1e-6, 1e-3]:
        x0_candidates.append([
            mu0,
            ah.inv_softplus(sigma0),
            math.log(q0_val / (1.0 - q0_val))
        ])

    # 2. Paper_V1 default
    m_hat = float(series_h.sum() / series_n.sum())
    mu_crude = float(norm.ppf(np.clip(m_hat, 1e-6, 1.0 - 1e-6)))
    x0_candidates.append([
        mu_crude,
        ah.inv_softplus(0.5),
        math.log(0.01 / 0.99)
    ])

    best_res = None
    best_nll = float("inf")

    for x0 in x0_candidates:
        try:
            res = minimize(
                ah.neg_ll_hier_LD,
                x0=np.array(x0, dtype=float),
                args=(n_arr, h_arr, c_arr, gh_x, gh_w),
                method="Nelder-Mead",
                options={"maxiter": maxiter}
            )
            if res.fun < best_nll:
                best_nll = float(res.fun)
                best_res = res
        except Exception:
            continue

    mu_hat = float(best_res.x[0])
    sigma_hat = float(ah.softplus(best_res.x[1]))
    q_hat = float(ah.sigmoid(best_res.x[2]))

    return {
        "mu": mu_hat,
        "sigma": sigma_hat,
        "q": q_hat,
        "nll": best_nll,
        "success": bool(best_res.success),
        "nit": int(best_res.nit),
        "message": str(best_res.message),
        "GH_N": int(GH_N),
    }

def nll_probit_normal_mixture_physical(mu1, mu2, sigma1, sigma2, pi, n_arr, h_arr, c_arr, GH_N=100):
    """
    Direct evaluation of Negative Log-Likelihood for 2-component Gaussian mixture environmental NULL (M3)
    using physical parameters (mu1, mu2, sigma1, sigma2, pi).
    y_t ~ pi N(mu1, sigma1^2) + (1-pi) N(mu2, sigma2^2)
    p_t = Phi(y_t)
    L_t | p_t, n_t ~ Binomial(n_t, p_t)
    """
    if not (0.0 < pi < 1.0) or sigma1 <= 0.0 or sigma2 <= 0.0:
        return float("inf")

    x, w = hermgauss(int(GH_N))
    w_nodes = w / math.sqrt(math.pi)
    w_nodes = w_nodes / np.sum(w_nodes)
    log_w = np.log(w_nodes)

    # Component 1
    y1 = mu1 + math.sqrt(2.0) * sigma1 * x
    p1 = np.clip(norm.cdf(y1), 1e-15, 1.0 - 1e-15)
    log_p1 = np.log(p1)
    log_1mp1 = np.log1p(-p1)

    # Component 2
    y2 = mu2 + math.sqrt(2.0) * sigma2 * x
    p2 = np.clip(norm.cdf(y2), 1e-15, 1.0 - 1e-15)
    log_p2 = np.log(p2)
    log_1mp2 = np.log1p(-p2)

    log_pi = math.log(pi)
    log_1mpi = math.log(1.0 - pi)

    tot_nll = 0.0
    for n, h, cnt in zip(n_arr, h_arr, c_arr):
        log_c = log_comb(int(n), int(h))

        terms1 = log_c + h * log_p1 + (n - h) * log_1mp1
        log_I1 = scalar_or_array_logsumexp(log_w + terms1)

        terms2 = log_c + h * log_p2 + (n - h) * log_1mp2
        log_I2 = scalar_or_array_logsumexp(log_w + terms2)

        log_prob = logsumexp_two_scalars(log_pi + log_I1, log_1mpi + log_I2)
        tot_nll -= cnt * log_prob

    return tot_nll

def test_collapsed_mixture_equivalence(series_n, series_h, GH_N=100):
    """
    Test direct collapsed mixture boundary (mu1 = mu2 = mu_M0, sigma1 = sigma2 = sigma_M0, pi = 0.5)
    against M0 NLL without using optimizer transformation.
    """
    n_arr, h_arr, c_arr = ah.make_grouped_nh(series_n, series_h)
    fit_m0 = fit_probit_normal_null(n_arr, h_arr, c_arr, GH_N=GH_N)
    mu0, sigma0 = fit_m0["mu"], fit_m0["sigma"]
    nll_m0 = fit_m0["nll"]

    nll_collapsed = nll_probit_normal_mixture_physical(
        mu0, mu0, sigma0, sigma0, 0.5, n_arr, h_arr, c_arr, GH_N=GH_N
    )
    delta = abs(nll_collapsed - nll_m0)
    return {
        "mu0": mu0,
        "sigma0": sigma0,
        "nll_m0": nll_m0,
        "nll_collapsed": nll_collapsed,
        "delta": delta,
        "collapsed_ok": bool(delta < 1e-6)
    }

def unconstrained_to_mixture_params(params):
    """
    Convert unconstrained vector [mu1, theta_dmu, theta_sig1, theta_sig2, theta_pi]
    to physical parameters (mu1, mu2, sigma1, sigma2, pi) with mu1 < mu2.
    """
    mu1 = float(params[0])
    dmu = float(ah.softplus(params[1]))
    mu2 = mu1 + dmu
    sigma1 = float(ah.softplus(params[2]))
    sigma2 = float(ah.softplus(params[3]))
    pi = float(ah.sigmoid(params[4]))
    return mu1, mu2, sigma1, sigma2, pi

def nll_probit_normal_mixture_unconstrained(params, n_arr, h_arr, c_arr, GH_N=100):
    mu1, mu2, sigma1, sigma2, pi = unconstrained_to_mixture_params(params)
    return nll_probit_normal_mixture_physical(mu1, mu2, sigma1, sigma2, pi, n_arr, h_arr, c_arr, GH_N=GH_N)

def fit_probit_normal_mixture_multistart(series_n, series_h, dataset_class="ALL", GH_N=100, maxiter=6000):
    """
    Multi-start optimization for M3 interior parameters + comparison against collapsed M0 boundary.
    Saves full per-start diagnostics to list of dicts.
    """
    n_arr, h_arr, c_arr = ah.make_grouped_nh(series_n, series_h)

    # 1. Collapsed M0 boundary evaluation
    collapsed_audit = test_collapsed_mixture_equivalence(series_n, series_h, GH_N=GH_N)
    mu0, sigma0 = collapsed_audit["mu0"], collapsed_audit["sigma0"]
    nll_m0 = collapsed_audit["nll_m0"]

    # 2. Build multi-start candidate physical initial points (interior mu1 < mu2)
    start_grid = []

    # Symmetric / asymmetric pi values
    pi_list = [0.1, 0.3, 0.5, 0.7, 0.9]
    # Separation around mu0
    mu_offsets = [0.1, 0.5, 1.0, 1.5]

    for pi_init in pi_list:
        for offset in mu_offsets:
            mu1_init = mu0 - offset / 2.0
            mu2_init = mu0 + offset / 2.0
            for sig_factor1 in [0.7, 1.0, 1.3]:
                for sig_factor2 in [0.7, 1.0, 1.3]:
                    sig1_init = max(1e-3, sigma0 * sig_factor1)
                    sig2_init = max(1e-3, sigma0 * sig_factor2)
                    start_grid.append((pi_init, mu1_init, sig1_init, mu2_init, sig2_init))

    diagnostics = []
    best_interior_res = None
    best_interior_nll = float("inf")
    best_interior_params = None

    for start_id, (pi_i, mu1_i, sig1_i, mu2_i, sig2_i) in enumerate(start_grid, start=1):
        x0 = [
            mu1_i,
            ah.inv_softplus(mu2_i - mu1_i),
            ah.inv_softplus(sig1_i),
            ah.inv_softplus(sig2_i),
            math.log(pi_i / (1.0 - pi_i))
        ]

        try:
            res = minimize(
                nll_probit_normal_mixture_unconstrained,
                x0=np.array(x0, dtype=float),
                args=(n_arr, h_arr, c_arr, GH_N),
                method="Nelder-Mead",
                options={"maxiter": maxiter}
            )
            final_mu1, final_mu2, final_sig1, final_sig2, final_pi = unconstrained_to_mixture_params(res.x)
            nll_val = float(res.fun)
            succ = bool(res.success)
            nit_val = int(res.nit)
        except Exception:
            final_mu1, final_mu2, final_sig1, final_sig2, final_pi = mu1_i, mu2_i, sig1_i, sig2_i, pi_i
            nll_val = float("inf")
            succ = False
            nit_val = 0

        diag_entry = {
            "class": dataset_class,
            "start_id": start_id,
            "initial_pi": pi_i,
            "initial_mu1": mu1_i,
            "initial_sigma1": sig1_i,
            "initial_mu2": mu2_i,
            "initial_sigma2": sig2_i,
            "final_pi": final_pi,
            "final_mu1": final_mu1,
            "final_sigma1": final_sig1,
            "final_mu2": final_mu2,
            "final_sigma2": final_sig2,
            "nll": nll_val,
            "success": succ,
            "nit": nit_val
        }
        diagnostics.append(diag_entry)

        if nll_val < best_interior_nll:
            best_interior_nll = nll_val
            best_interior_params = (final_mu1, final_mu2, final_sig1, final_sig2, final_pi)
            best_interior_res = res

    # 3. Overall M3 best-found selection (comparison with collapsed boundary)
    if best_interior_nll < nll_m0 - 1e-6:
        best_nll = best_interior_nll
        mu1_best, mu2_best, sig1_best, sig2_best, pi_best = best_interior_params
        is_collapsed = False
    else:
        best_nll = nll_m0
        mu1_best, mu2_best = mu0, mu0
        sig1_best, sig2_best = sigma0, sigma0
        pi_best = 0.5
        is_collapsed = True

    return {
        "class": dataset_class,
        "nll_m0": nll_m0,
        "nll_best_interior": best_interior_nll,
        "nll": best_nll,
        "mu1": mu1_best,
        "mu2": mu2_best,
        "sigma1": sig1_best,
        "sigma2": sig2_best,
        "pi": pi_best,
        "is_collapsed_boundary": is_collapsed,
        "diagnostics": diagnostics
    }

