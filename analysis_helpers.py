"""Analysis helpers for summary statistics, hierarchical fits, and variance decomposition."""

import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from numpy.polynomial.hermite import hermgauss


CLASSES = ["ALL", "SG", "IG"]
GH_N_HIER = 100
MAXITER_IID = 8000
MAXITER_HIER = 6000
ALPHAS = [0.95, 0.99]
from models import (
    log_pmf_Torri,
    log_pmf_LD_1pt,
    log_pmf_Vas_1pt,
    normalize_logpmf,
    survival_from_pmf,
    tail_metrics_from_pmf,
    m_rho_from_pmf,
    moments_from_pmf,
)


def summarize_nt_default(df, label):
    n = df[label].to_numpy(dtype=float)
    L = df["D_" + label].to_numpy(dtype=float)

    return {
        "Class": label,
        "Mean n": np.mean(n),
        "Mean L": np.mean(L),
        "Mean default rate": np.mean(L / n),
        "Total default rate": np.sum(L) / np.sum(n),
    }




def summarize_nt_default_by_period(df, label):
    out = []
    for period in ["1950-1979", "1980-2023"]:
        g = df[df["period"] == period]

        n = g[label].to_numpy(dtype=float)
        L = g["D_" + label].to_numpy(dtype=float)

        out.append({
            "Class": label,
            "Period": period,
            "Mean n": np.mean(n),
            "Mean L": np.mean(L),
            "Mean default rate": np.mean(L / n),
            "Total default rate": np.sum(L) / np.sum(n),
        })
    return out




def fit_one_Torri(series_n, series_h, x0=None, maxiter=5000):
    data = np.c_[series_n.to_numpy(dtype=int), series_h.to_numpy(dtype=int)]

    def neg_ll_torri(theta):
        # unconstrained -> (0,1) via sigmoid
        p = 1.0/(1.0 + np.exp(-theta[0]))
        u = 1.0/(1.0 + np.exp(-theta[1]))
        v = 1.0/(1.0 + np.exp(-theta[2]))

        ll = 0.0
        for n, h in data:
            lp = log_pmf_Torri(int(n), int(h), float(p), float(u), float(v))
            if not np.isfinite(lp):
                return 1e100
            ll += lp
        return -ll

    # initial guess (logit)
    if x0 is None:
        x0 = np.array([
            math.log(0.01/0.99),  # p ~ 1%
            math.log(0.9/0.1),    # u ~ 0.9
            math.log(0.1/0.9),    # v ~ 0.1
        ], dtype=float)

    res = minimize(neg_ll_torri, x0, method="Nelder-Mead", options={"maxiter": maxiter})

    th = res.x
    p = 1.0/(1.0 + np.exp(-th[0]))
    u = 1.0/(1.0 + np.exp(-th[1]))
    v = 1.0/(1.0 + np.exp(-th[2]))
    return {
        "p": float(p), "u": float(u), "v": float(v),
        "nll": float(res.fun),
        "success": bool(res.success),
        "nit": int(res.nit),
        "message": str(res.message),
    }




def fit_one_LD(series_n, series_h):
    tmp = pd.DataFrame({"n": series_n.to_numpy(dtype=int),
                        "h": series_h.to_numpy(dtype=int)})
    g = tmp.value_counts().reset_index(name="cnt")
    n_arr = g["n"].to_numpy(dtype=int)
    h_arr = g["h"].to_numpy(dtype=int)
    c_arr = g["cnt"].to_numpy(dtype=int)

    def neg_ll_ld(theta):
        p = 1/(1+np.exp(-theta[0]))
        q = 1/(1+np.exp(-theta[1]))

        ll = 0.0
        for n, h, c in zip(n_arr, h_arr, c_arr):
            lp = log_pmf_LD_1pt(int(n), int(h), float(p), float(q))
            if not np.isfinite(lp):
                return 1e100
            ll += c * lp
        return -ll

    x0 = np.array([math.log(0.01/0.99), math.log(0.02/0.98)])
    res = minimize(neg_ll_ld, x0, method="Nelder-Mead", options={"maxiter":8000})

    th = res.x
    p = 1/(1+np.exp(-th[0]))
    q = 1/(1+np.exp(-th[1]))
    return {
        "p": float(p), "q": float(q),
        "nll": float(res.fun),
        "success": bool(res.success),
        "nit": int(res.nit),
        "message": str(res.message),
    }




def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))




def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)




def inv_softplus(s):
    s = float(s)
    if s <= 0:
        raise ValueError("softplus inverse requires s > 0")
    return s + np.log(-np.expm1(-s))




def _clean_nh(series_n, series_h):
    tmp = pd.DataFrame({
        "n": series_n.to_numpy(dtype=int),
        "h": series_h.to_numpy(dtype=int),
    })
    tmp = tmp[(tmp["n"] > 0) & (tmp["h"] >= 0) & (tmp["h"] <= tmp["n"])].copy()
    return tmp




def make_grouped_nh(series_n, series_h):
    tmp = _clean_nh(series_n, series_h)
    g = tmp.value_counts().reset_index(name="cnt")
    n_arr = g["n"].to_numpy(dtype=int)
    h_arr = g["h"].to_numpy(dtype=int)
    c_arr = g["cnt"].to_numpy(dtype=int)
    return n_arr, h_arr, c_arr




def _clean_nL(df, key):
    n_vals = df[key].astype(int).to_numpy()
    L_vals = df["D_" + key].astype(int).to_numpy()
    mask = (
        (n_vals > 0) &
        (L_vals >= 0) &
        (L_vals <= n_vals) &
        np.isfinite(n_vals) &
        np.isfinite(L_vals)
    )
    return n_vals[mask].astype(int), L_vals[mask].astype(int)




def scaled_counts(n_vals, L_vals, nbar):
    n_vals = np.asarray(n_vals, dtype=float)
    L_vals = np.asarray(L_vals, dtype=float)
    Ls = np.rint(L_vals * (float(nbar) / n_vals)).astype(int)
    Ls = np.maximum(Ls, 0)
    return Ls




def gh_normal_weights(mu, sigma, GH_N):
    x, w = hermgauss(int(GH_N))
    y = mu + math.sqrt(2.0) * sigma * x
    wn = w / math.sqrt(math.pi)
    wn = wn / np.sum(wn)
    return y, wn




def aggregate_m_q_to_rho(m, q):
    den = m * (1.0 - m)
    if den <= 1e-15:
        return float("nan")
    return (q - m * m) / den




def torri_m_q_for_n(n, p, u, v):
    b = (1.0 - p) * (1.0 - u)
    s = 1.0 - p * v
    a = p * (1.0 - v)
    m = p + b * (1.0 - s**(n - 1)) if n >= 1 else 0.0
    q_pair = (p + b)**2 - b * (2.0 * a + b) * (s**(n - 2)) if n >= 2 else m * m
    return float(m), float(q_pair)




def build_compare_decomp_table(out_var_hier_ld, out_var_hier_torri):
    cols = [
        "class",
        "Share_iid_within_model",
        "Share_infect_within_model",
        "Share_pt_within_model",
        "Iid_over_VarData",
        "Infect_over_VarData",
        "Pt_over_VarData",
    ]
    df_ld = out_var_hier_ld[cols].copy()
    df_ld["model"] = "Lo–Davis"
    df_t = out_var_hier_torri[cols].copy()
    df_t["model"] = "Torri"
    df_cmp = pd.concat([df_ld, df_t], axis=0, ignore_index=True)
    return df_cmp.sort_values(["class", "model"]).reset_index(drop=True)




def empirical_pmf_from_counts(Ls):
    Ls = np.asarray(Ls, dtype=int)
    hmax = int(Ls.max())
    h = np.arange(hmax + 1)
    hist, _ = np.histogram(Ls, bins=np.arange(hmax + 2) - 0.5, density=True)
    return h, hist




def model_pmf_fixed_n(model, key, nbar, hmax):
    h = np.arange(hmax + 1, dtype=int)

    if model == "Torri":
        p,u,v = fits_torri[key]["p"], fits_torri[key]["u"], fits_torri[key]["v"]
        lp = np.full(hmax + 1, -np.inf, dtype=float)
        h_ok = h[h <= nbar]
        lp[h_ok] = np.array([log_pmf_Torri(int(nbar), int(x), float(p), float(u), float(v)) for x in h_ok], dtype=float)
        pmf = normalize_logpmf(lp)
        return pmf

    if model == "LD":
        p,q = fits_ld[key]["p"], fits_ld[key]["q"]
        lp = np.full(hmax + 1, -np.inf, dtype=float)
        h_ok = h[h <= nbar]
        lp[h_ok] = np.array([log_pmf_LD_1pt(int(nbar), int(x), float(p), float(q)) for x in h_ok], dtype=float)
        pmf = normalize_logpmf(lp)
        return pmf

    if model == "Vas":
        p,rhoA = fits_vas[key]["p"], fits_vas[key]["rho_A"]
        lp = np.full(hmax + 1, -np.inf, dtype=float)
        h_ok = h[h <= nbar]
        lp[h_ok] = np.array([log_pmf_Vas_1pt(int(nbar), int(x), float(p), float(rhoA)) for x in h_ok], dtype=float)
        pmf = normalize_logpmf(lp)
        return pmf




def data_m_q_hat(n_vals, L_vals):
    """Pair-weighted estimates from annual (n_t, L_t)."""
    n = np.asarray(n_vals, dtype=float)
    L = np.asarray(L_vals, dtype=float)

    # m-hat: weighted by n_t (natural for mean default rate over obligors)
    sum_n = np.sum(n)
    m_hat = np.sum(L) / sum_n

    # q-hat: weighted by n_t(n_t-1) using E[L(L-1)] identity
    w_pair = n * (n - 1.0)
    sum_w = np.sum(w_pair)
    if sum_w <= 0:
        q_hat = float("nan")
    else:
        q_hat = np.sum(L * (L - 1.0)) / sum_w

    rho_hat = aggregate_m_q_to_rho(m_hat, q_hat)
    return float(m_hat), float(q_hat), float(rho_hat)




def model_aggregate_m_q(n_vals, m_q_fn):
    """
    Aggregate model-implied m(n), q(n) over observed n_t with symmetric weights.
    m_q_fn(n) -> (m_n, q_n) where
      m_n = E[Z_i] (default prob)
      q_n = P(Z_i=1, Z_j=1) for i!=j
    """
    n = np.asarray(n_vals, dtype=int)
    uniq, cnt = np.unique(n, return_counts=True)

    # weights for m: sum over years of n_t
    sum_n = float(np.sum(uniq * cnt))

    # weights for q: sum over years of n_t(n_t-1)
    sum_pair = float(np.sum(uniq * (uniq - 1) * cnt))

    num_m = 0.0
    num_q = 0.0

    for nn, cc in zip(uniq, cnt):
        if nn <= 0:
            continue
        m_n, q_n = m_q_fn(int(nn))
        num_m += float(cc) * float(nn) * float(m_n)
        if nn >= 2:
            num_q += float(cc) * float(nn) * float(nn - 1) * float(q_n)

    m = num_m / sum_n if sum_n > 0 else float("nan")
    q = num_q / sum_pair if sum_pair > 0 else float("nan")
    rho = aggregate_m_q_to_rho(m, q)
    return float(m), float(q), float(rho)




def torri_m_q_fn(p, u, v):
    # Use your existing closed-form moments_from_params_Torri(n,p,u,v) -> (m, rho)
    # Then convert to q = m^2 + rho*m(1-m)
    def _fn(n):
        m_n, rho_n = moments_from_params_Torri(n, p, u, v)
        q_n = m_n*m_n + rho_n*m_n*(1.0-m_n) if np.isfinite(rho_n) else float("nan")
        return float(m_n), float(q_n)
    return _fn




def ld_m_q_fn(p, q):
    # Use your existing closed-form m_rho_LD(n,p,q) -> (m, rho)
    def _fn(n):
        m_n, rho_n = m_rho_LD(n, p, q)
        q_n = m_n*m_n + rho_n*m_n*(1.0-m_n) if np.isfinite(rho_n) else float("nan")
        return float(m_n), float(q_n)
    return _fn




def vas_m_q_fn(p, rhoA):
    # Vasicek: m is constant = p.
    # Pair default prob q = Phi_2(a,a; rhoA) where a = Phi^{-1}(p)
    a = norm.ppf(p)
    q_pair = phi2_aa(a, rhoA)
    def _fn(n):
        return float(p), float(q_pair)
    return _fn




def ld_m_q_for_n(n, p, q_param):
    # m = 1 - (1-p)(1-qp)^(n-1)
    m = 1.0 - (1.0 - p) * (1.0 - q_param * p)**(n - 1) if n >= 1 else 0.0
    if n < 2:
        return float(m), float(m*m)

    # cov formula (from your m_rho_LD docstring)
    termA = (1.0 - 2.0*p*q_param + p*(q_param**2))**(n - 2)
    termB = (1.0 - p*q_param)**(2*(n - 1))
    cov = (1.0 - p)**2 * (termA - termB)

    # q = m^2 + cov
    q_pair = m*m + cov
    return float(m), float(q_pair)




def vas_m_q_for_n(n, p, rhoA):
    m = float(p)
    a = norm.ppf(m)
    q_pair = float(phi2_aa(a, float(rhoA)))  # P(X<=a, Y<=a)
    return float(m), float(q_pair)




def variance_explained_ratio_3models(
    df,
    fits_torri,   # dict: {"ALL":{"p":..,"u":..,"v":..}, ...}
    fits_ld,      # dict: {"ALL":{"p":..,"q":..}, ...}
    fits_vas,     # dict: {"ALL":{"p":..,"rho_A":..,"GH_N":..}, ...}  (GH_N optional)
    keys=("ALL","SG","IG")
):
    rows = []

    for key in keys:
        n, L = _clean_nL(df, key)
        if len(L) < 2:
            continue

        # DATA: annual-count variance across years
        var_data = float(np.var(L, ddof=1))

        # DATA: obligor-year weighted mean default rate
        m_hat = float(L.sum() / n.sum())

        # crude iid baseline (uses global m_hat)
        var_iid_data = float(np.mean(n * m_hat * (1 - m_hat)))
        var_corr_data = var_data - var_iid_data

        # iterate models
        for model in ["Torri", "Lo–Davis", "Vasicek"]:
            var_iid_list = []
            var_corr_list = []

            if model == "Torri":
                p = float(fits_torri[key]["p"])
                u = float(fits_torri[key]["u"])
                v = float(fits_torri[key]["v"])
                for nn in n:
                    m_nn, q_nn = torri_m_q_for_n(int(nn), p, u, v)
                    cov_nn = q_nn - m_nn**2
                    var_iid_list.append(nn * m_nn * (1 - m_nn))
                    var_corr_list.append(nn * (nn - 1) * cov_nn)

            elif model == "Lo–Davis":
                p = float(fits_ld[key]["p"])
                q_param = float(fits_ld[key]["q"])
                for nn in n:
                    m_nn, q_nn = ld_m_q_for_n(int(nn), p, q_param)
                    cov_nn = q_nn - m_nn**2
                    var_iid_list.append(nn * m_nn * (1 - m_nn))
                    var_corr_list.append(nn * (nn - 1) * cov_nn)

            else:  # Vasicek
                p = float(fits_vas[key]["p"])
                rhoA = float(fits_vas[key]["rho_A"])
                for nn in n:
                    m_nn, q_nn = vas_m_q_for_n(int(nn), p, rhoA)
                    cov_nn = q_nn - m_nn**2
                    var_iid_list.append(nn * m_nn * (1 - m_nn))
                    var_corr_list.append(nn * (nn - 1) * cov_nn)

            var_iid_model = float(np.mean(var_iid_list))
            var_corr_model = float(np.mean(var_corr_list))
            var_model = var_iid_model + var_corr_model

            # (B-like): corr contribution relative to total annual variance
            explained_vs_data = (var_corr_model / var_data) if var_data > 0 else np.nan

            # (A): corr share within conditional variance
            corr_share_cond = (var_corr_model / var_model) if var_model > 0 else np.nan

            rows.append({
                "class": key,
                "model": model,
                "Var_data": var_data,
                "m_hat_data": m_hat,
                "Var_iid_data": var_iid_data,
                "Var_corr_data": var_corr_data,
                "Var_model": var_model,
                "Var_iid_model": var_iid_model,
                "Var_corr_model": var_corr_model,
                "Explained_corr_over_VarData": explained_vs_data,     # can exceed 1
                "Corr_share_within_condVar": corr_share_cond,         # in [0,1] ideally
            })

    return pd.DataFrame(rows)




def log_pmf_hier_LD_1pt(n, h, mu, sigma, q, gh_x, gh_w):
    """
    log P(L=h | n, mu, sigma, q)
      = log ∫ P_LD(L=h | n, p=Phi(y), q) * N(y | mu, sigma^2) dy
      = log [ (1/sqrt(pi)) sum_j w_j * P_LD(h | n, Phi(mu + sqrt(2)sigma x_j), q) ]
    """
    if not (0 <= h <= n):
        return -np.inf
    if sigma <= 0:
        return -np.inf
    if not (0.0 < q < 1.0):
        return -np.inf

    # GH approximation
    log_terms = np.full(len(gh_w), -np.inf, dtype=float)

    for j, (xj, wj) in enumerate(zip(gh_x, gh_w)):
        y = mu + math.sqrt(2.0) * sigma * float(xj)
        p = float(norm.cdf(y))

        # protect against exact 0/1
        p = float(np.clip(p, 1e-15, 1.0 - 1e-15))

        lp = log_pmf_LD_1pt(int(n), int(h), p, float(q))
        if np.isfinite(lp):
            log_terms[j] = math.log(float(wj)) + lp

    out = logsumexp_vec(log_terms) - 0.5 * math.log(math.pi)
    return float(out)




def neg_ll_hier_LD(theta, n_arr, h_arr, c_arr, gh_x, gh_w):
    """
    theta = (mu_raw, sigma_raw, q_raw)
      mu    = mu_raw
      sigma = softplus(sigma_raw)
      q     = sigmoid(q_raw)
    """
    mu = float(theta[0])
    sigma = float(softplus(theta[1]))
    q = float(sigmoid(theta[2]))

    ll = 0.0
    for n, h, c in zip(n_arr, h_arr, c_arr):
        lp = log_pmf_hier_LD_1pt(int(n), int(h), mu, sigma, q, gh_x, gh_w)
        if not np.isfinite(lp):
            return 1e100
        ll += int(c) * lp

    return float(-ll)




def fit_one_hier_LD(
    series_n,
    series_h,
    GH_N=20,
    x0=None,
    maxiter=5000,
    method="Nelder-Mead",
):
    """
    Hierarchical LD MLE.

    Parameters
    ----------
    series_n, series_h : pandas Series
        Observed n_t and L_t
    GH_N : int
        Number of Gauss-Hermite nodes
    x0 : array-like or None
        Initial values for unconstrained params (mu_raw, sigma_raw, q_raw)
    """
    n_arr, h_arr, c_arr = make_grouped_nh(series_n, series_h)
    gh_x, gh_w = hermgauss(GH_N)

    # crude initial values from overall default rate
    m_hat = float(series_h.sum() / series_n.sum())
    mu0 = float(norm.ppf(np.clip(m_hat, 1e-6, 1.0 - 1e-6)))
    sigma0 = 0.5
    q0 = 0.01

    if x0 is None:
        x0 = np.array([
            mu0,
            inv_softplus(sigma0),
            math.log(q0 / (1.0 - q0)),
        ], dtype=float)

    res = minimize(
        neg_ll_hier_LD,
        x0=x0,
        args=(n_arr, h_arr, c_arr, gh_x, gh_w),
        method=method,
        options={"maxiter": maxiter},
    )

    mu_hat = float(res.x[0])
    sigma_hat = float(softplus(res.x[1]))
    q_hat = float(sigmoid(res.x[2]))

    return {
        "mu": mu_hat,
        "sigma": sigma_hat,
        "q": q_hat,
        "nll": float(res.fun),
        "success": bool(res.success),
        "nit": int(res.nit),
        "message": str(res.message),
        "GH_N": int(GH_N),
        "m_hat_init": m_hat,
    }




def posterior_mean_pt_hier_LD(n, h, mu, sigma, q, gh_x, gh_w):
    """
    E[p_t | L_t=h, n_t=n, mu, sigma, q] by GH approximation
    """
    logw = np.full(len(gh_w), -np.inf, dtype=float)
    pvals = np.zeros(len(gh_w), dtype=float)

    for j, (xj, wj) in enumerate(zip(gh_x, gh_w)):
        y = mu + math.sqrt(2.0) * sigma * float(xj)
        p = float(norm.cdf(y))
        p = float(np.clip(p, 1e-15, 1.0 - 1e-15))
        lp = log_pmf_LD_1pt(int(n), int(h), p, float(q))
        if np.isfinite(lp):
            logw[j] = math.log(float(wj)) + lp
            pvals[j] = p

    z = logsumexp_vec(logw)
    if not np.isfinite(z):
        return np.nan

    w_norm = np.exp(logw - z)
    return float(np.sum(w_norm * pvals))




def add_posterior_mean_pt_column(df_in, key, fit_out):
    """
    Add posterior mean p_t column for one class.
    """
    gh_x, gh_w = hermgauss(int(fit_out["GH_N"]))
    mu = float(fit_out["mu"])
    sigma = float(fit_out["sigma"])
    q = float(fit_out["q"])

    out = df_in.copy()
    p_post = []
    for n, h in zip(out[key].astype(int), out["D_" + key].astype(int)):
        p_post.append(posterior_mean_pt_hier_LD(int(n), int(h), mu, sigma, q, gh_x, gh_w))
    out["p_post_mean_" + key] = p_post
    return out




def variance_decomp_hier_ld_fixed_nbar(
    df,
    NBAR,
    fits_hier_ld,
    keys=("ALL", "SG", "IG"),
):
    rows = []

    for key in keys:
        nbar = int(NBAR[key])
        n_vals, L_vals = _clean_nL(df, key)
        L_scaled = scaled_counts(n_vals, L_vals, nbar)

        # data-side quantities under the same fixed-nbar scaling
        var_data_scaled = float(np.var(L_scaled, ddof=1)) if len(L_scaled) >= 2 else np.nan
        m_hat_data = float(np.sum(L_vals) / np.sum(n_vals))
        var_iid_data_scaled = float(nbar * m_hat_data * (1.0 - m_hat_data))
        var_corr_data_scaled = var_data_scaled - var_iid_data_scaled if np.isfinite(var_data_scaled) else np.nan

        # fitted hyperparameters
        mu = float(fits_hier_ld[key]["mu"])
        sigma = float(fits_hier_ld[key]["sigma"])
        q_param = float(fits_hier_ld[key]["q"])
        GH_N = int(fits_hier_ld[key]["GH_N"])

        # GH quadrature for p_t = Phi(y_t), y_t ~ N(mu, sigma^2)
        y_nodes, w_nodes = gh_normal_weights(mu, sigma, GH_N)
        p_nodes = norm.cdf(y_nodes)
        p_nodes = np.clip(p_nodes, 1e-15, 1.0 - 1e-15)

        # model-implied m_t, q_t at fixed nbar for each quadrature node
        m_nodes = np.empty_like(p_nodes, dtype=float)
        q_nodes = np.empty_like(p_nodes, dtype=float)

        for j, p_t in enumerate(p_nodes):
            m_t, q_t = ld_m_q_for_n(nbar, float(p_t), q_param)
            m_nodes[j] = m_t
            q_nodes[j] = q_t

        # decomposition
        var_iid_model = float(np.sum(w_nodes * (nbar * m_nodes * (1.0 - m_nodes))))
        var_infect_model = float(np.sum(w_nodes * (nbar * (nbar - 1) * (q_nodes - m_nodes**2))))
        var_pt_model = float(np.sum(w_nodes * (nbar * m_nodes)**2) - (np.sum(w_nodes * (nbar * m_nodes)))**2)

        var_model_total = var_iid_model + var_infect_model + var_pt_model

        # additional summaries
        m_bar = float(np.sum(w_nodes * m_nodes))
        q_bar = float(np.sum(w_nodes * q_nodes))
        rho_bar = (q_bar - m_bar*m_bar) / (m_bar * (1.0 - m_bar)) if m_bar * (1.0 - m_bar) > 1e-15 else np.nan

        rows.append({
            "class": key,
            "nbar": nbar,
            "mu": mu,
            "sigma": sigma,
            "q_param": q_param,
            "GH_N": GH_N,
            "Var_data_scaled": var_data_scaled,
            "m_hat_data": m_hat_data,
            "Var_iid_data_scaled": var_iid_data_scaled,
            "Var_corr_data_scaled": var_corr_data_scaled,
            "m_bar_model": m_bar,
            "q_bar_model": q_bar,
            "rho_bar_model": rho_bar,
            "Var_iid_model": var_iid_model,
            "Var_infect_model": var_infect_model,
            "Var_pt_model": var_pt_model,
            "Var_model_total": var_model_total,
            "Share_iid_within_model": (var_iid_model / var_model_total) if var_model_total > 0 else np.nan,
            "Share_infect_within_model": (var_infect_model / var_model_total) if var_model_total > 0 else np.nan,
            "Share_pt_within_model": (var_pt_model / var_model_total) if var_model_total > 0 else np.nan,
            "Infect_over_VarData": (var_infect_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
            "Pt_over_VarData": (var_pt_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
            "Iid_over_VarData": (var_iid_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
        })

    return pd.DataFrame(rows)




def log_pmf_hier_Torri_1pt(n, h, mu, sigma, u, v, gh_x, gh_w):
    """
    log P(L=h | n, mu, sigma, u, v)
      = log ∫ P_Torri(L=h | n, p=Phi(y), u, v) * N(y | mu, sigma^2) dy
      = log [ (1/sqrt(pi)) sum_j w_j * P_Torri(h | n, Phi(mu + sqrt(2)sigma x_j), u, v) ]
    """
    if not (0 <= h <= n):
        return -np.inf
    if sigma <= 0:
        return -np.inf
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        return -np.inf

    log_terms = np.full(len(gh_w), -np.inf, dtype=float)

    for j, (xj, wj) in enumerate(zip(gh_x, gh_w)):
        y = mu + math.sqrt(2.0) * sigma * float(xj)
        p = float(norm.cdf(y))
        p = float(np.clip(p, 1e-15, 1.0 - 1e-15))

        lp = log_pmf_Torri(int(n), int(h), p, float(u), float(v))
        if np.isfinite(lp):
            log_terms[j] = math.log(float(wj)) + lp

    out = logsumexp_vec(log_terms) - 0.5 * math.log(math.pi)
    return float(out)




def neg_ll_hier_Torri(theta, n_arr, h_arr, c_arr, gh_x, gh_w):
    """
    theta = (mu_raw, sigma_raw, u_raw, v_raw)
      mu    = mu_raw
      sigma = softplus(sigma_raw)
      u     = sigmoid(u_raw)
      v     = sigmoid(v_raw)
    """
    mu = float(theta[0])
    sigma = float(softplus(theta[1]))
    u = float(sigmoid(theta[2]))
    v = float(sigmoid(theta[3]))

    ll = 0.0
    for n, h, c in zip(n_arr, h_arr, c_arr):
        lp = log_pmf_hier_Torri_1pt(int(n), int(h), mu, sigma, u, v, gh_x, gh_w)
        if not np.isfinite(lp):
            return 1e100
        ll += int(c) * lp

    return float(-ll)




def fit_one_hier_Torri(
    series_n,
    series_h,
    GH_N=20,
    x0=None,
    maxiter=7000,
    method="Nelder-Mead",
):
    """
    Hierarchical Torri MLE.
    """
    n_arr, h_arr, c_arr = make_grouped_nh(series_n, series_h)
    gh_x, gh_w = hermgauss(GH_N)

    # crude initial values
    m_hat = float(series_h.sum() / series_n.sum())
    mu0 = float(norm.ppf(np.clip(m_hat, 1e-6, 1.0 - 1e-6)))
    sigma0 = 0.5
    u0 = 0.9
    v0 = 0.01

    if x0 is None:
        x0 = np.array([
            mu0,
            inv_softplus(sigma0),
            math.log(u0 / (1.0 - u0)),
            math.log(v0 / (1.0 - v0)),
        ], dtype=float)

    res = minimize(
        neg_ll_hier_Torri,
        x0=x0,
        args=(n_arr, h_arr, c_arr, gh_x, gh_w),
        method=method,
        options={"maxiter": maxiter},
    )

    mu_hat = float(res.x[0])
    sigma_hat = float(softplus(res.x[1]))
    u_hat = float(sigmoid(res.x[2]))
    v_hat = float(sigmoid(res.x[3]))

    return {
        "mu": mu_hat,
        "sigma": sigma_hat,
        "u": u_hat,
        "v": v_hat,
        "nll": float(res.fun),
        "success": bool(res.success),
        "nit": int(res.nit),
        "message": str(res.message),
        "GH_N": int(GH_N),
        "m_hat_init": m_hat,
    }




def posterior_mean_pt_hier_Torri(n, h, mu, sigma, u, v, gh_x, gh_w):
    logw = np.full(len(gh_w), -np.inf, dtype=float)
    pvals = np.zeros(len(gh_w), dtype=float)

    for j, (xj, wj) in enumerate(zip(gh_x, gh_w)):
        y = mu + math.sqrt(2.0) * sigma * float(xj)
        p = float(norm.cdf(y))
        p = float(np.clip(p, 1e-15, 1.0 - 1e-15))
        lp = log_pmf_Torri(int(n), int(h), p, float(u), float(v))
        if np.isfinite(lp):
            logw[j] = math.log(float(wj)) + lp
            pvals[j] = p

    z = logsumexp_vec(logw)
    if not np.isfinite(z):
        return np.nan

    w_norm = np.exp(logw - z)
    return float(np.sum(w_norm * pvals))




def add_posterior_mean_pt_column_torri(df_in, key, fit_out):
    gh_x, gh_w = hermgauss(int(fit_out["GH_N"]))
    mu = float(fit_out["mu"])
    sigma = float(fit_out["sigma"])
    u = float(fit_out["u"])
    v = float(fit_out["v"])

    out = df_in.copy()
    p_post = []
    for n, h in zip(out[key].astype(int), out["D_" + key].astype(int)):
        p_post.append(posterior_mean_pt_hier_Torri(int(n), int(h), mu, sigma, u, v, gh_x, gh_w))
    out["p_post_mean_" + key] = p_post
    return out




def variance_decomp_hier_torri_fixed_nbar(
    df,
    NBAR,
    fits_hier_torri,
    keys=("ALL", "SG", "IG"),
):
    rows = []

    for key in keys:
        nbar = int(NBAR[key])

        # data-side quantities under the same fixed-nbar scaling
        n_vals, L_vals = _clean_nL(df, key)
        L_scaled = scaled_counts(n_vals, L_vals, nbar)
        var_data_scaled = float(np.var(L_scaled, ddof=1)) if len(L_scaled) >= 2 else np.nan
        m_hat_data = float(np.sum(L_vals) / np.sum(n_vals))
        var_iid_data_scaled = float(nbar * m_hat_data * (1.0 - m_hat_data))
        var_corr_data_scaled = var_data_scaled - var_iid_data_scaled if np.isfinite(var_data_scaled) else np.nan

        # fitted hyperparameters
        mu = float(fits_hier_torri[key]["mu"])
        sigma = float(fits_hier_torri[key]["sigma"])
        u = float(fits_hier_torri[key]["u"])
        v = float(fits_hier_torri[key]["v"])
        GH_N = int(fits_hier_torri[key]["GH_N"])

        # GH quadrature for p_t = Phi(y_t), y_t ~ N(mu, sigma^2)
        y_nodes, w_nodes = gh_normal_weights(mu, sigma, GH_N)
        p_nodes = norm.cdf(y_nodes)
        p_nodes = np.clip(p_nodes, 1e-15, 1.0 - 1e-15)

        # model-implied m_t, q_t at fixed nbar for each quadrature node
        m_nodes = np.empty_like(p_nodes, dtype=float)
        q_nodes = np.empty_like(p_nodes, dtype=float)

        for j, p_t in enumerate(p_nodes):
            m_t, q_t = torri_m_q_for_n(nbar, float(p_t), u, v)
            m_nodes[j] = m_t
            q_nodes[j] = q_t

        # decomposition
        var_iid_model = float(np.sum(w_nodes * (nbar * m_nodes * (1.0 - m_nodes))))
        var_infect_model = float(np.sum(w_nodes * (nbar * (nbar - 1) * (q_nodes - m_nodes**2))))
        var_pt_model = float(np.sum(w_nodes * (nbar * m_nodes)**2) - (np.sum(w_nodes * (nbar * m_nodes)))**2)

        var_model_total = var_iid_model + var_infect_model + var_pt_model

        # additional summaries
        m_bar = float(np.sum(w_nodes * m_nodes))
        q_bar = float(np.sum(w_nodes * q_nodes))
        rho_bar = aggregate_m_q_to_rho(m_bar, q_bar)

        rows.append({
            "class": key,
            "nbar": nbar,
            "mu": mu,
            "sigma": sigma,
            "u": u,
            "v": v,
            "GH_N": GH_N,
            "Var_data_scaled": var_data_scaled,
            "m_hat_data": m_hat_data,
            "Var_iid_data_scaled": var_iid_data_scaled,
            "Var_corr_data_scaled": var_corr_data_scaled,
            "m_bar_model": m_bar,
            "q_bar_model": q_bar,
            "rho_bar_model": rho_bar,
            "Var_iid_model": var_iid_model,
            "Var_infect_model": var_infect_model,
            "Var_pt_model": var_pt_model,
            "Var_model_total": var_model_total,
            "Share_iid_within_model": (var_iid_model / var_model_total) if var_model_total > 0 else np.nan,
            "Share_infect_within_model": (var_infect_model / var_model_total) if var_model_total > 0 else np.nan,
            "Share_pt_within_model": (var_pt_model / var_model_total) if var_model_total > 0 else np.nan,
            "Iid_over_VarData": (var_iid_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
            "Infect_over_VarData": (var_infect_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
            "Pt_over_VarData": (var_pt_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
        })

    return pd.DataFrame(rows)




def pmf_vas_fixed_n(nbar, p, rho_A):
    h = np.arange(nbar + 1, dtype=int)
    lp = np.array([log_pmf_Vas_1pt(int(nbar), int(x), float(p), float(rho_A)) for x in h], dtype=float)
    return normalize_logpmf(lp)




def pmf_hier_ld_fixed_n(nbar, mu, sigma, q_param, GH_N):
    y_nodes, w_nodes = gh_normal_weights(mu, sigma, GH_N)
    p_nodes = np.clip(norm.cdf(y_nodes), 1e-15, 1.0 - 1e-15)

    h = np.arange(nbar + 1, dtype=int)
    mix = np.zeros(nbar + 1, dtype=float)

    for pj, wj in zip(p_nodes, w_nodes):
        lp = np.array([log_pmf_LD_1pt(int(nbar), int(x), float(pj), float(q_param)) for x in h], dtype=float)
        pmf = normalize_logpmf(lp)
        mix += wj * pmf

    mix = np.maximum(mix, 0.0)
    mix /= mix.sum()
    return mix




def pmf_hier_torri_fixed_n(nbar, mu, sigma, u, v, GH_N):
    y_nodes, w_nodes = gh_normal_weights(mu, sigma, GH_N)
    p_nodes = np.clip(norm.cdf(y_nodes), 1e-15, 1.0 - 1e-15)

    h = np.arange(nbar + 1, dtype=int)
    mix = np.zeros(nbar + 1, dtype=float)

    for pj, wj in zip(p_nodes, w_nodes):
        lp = np.array([log_pmf_Torri(int(nbar), int(x), float(pj), float(u), float(v)) for x in h], dtype=float)
        pmf = normalize_logpmf(lp)
        mix += wj * pmf

    mix = np.maximum(mix, 0.0)
    mix /= mix.sum()
    return mix




def model_pmf_fixed_n_hierarchical(model, key, nbar, hmax):
    if model == "LD":
        fit = fits_hier_ld[key]
        mu = float(fit["mu"])
        sigma = float(fit["sigma"])
        q_param = float(fit["q"])
        GH_N = int(fit["GH_N"])

    elif model == "Torri":
        fit = fits_hier_torri[key]
        mu = float(fit["mu"])
        sigma = float(fit["sigma"])
        u = float(fit["u"])
        v = float(fit["v"])
        GH_N = int(fit["GH_N"])

    else:
        raise ValueError("model must be 'LD' or 'Torri'")

    x, w = hermgauss(GH_N)
    y_nodes = mu + math.sqrt(2.0) * sigma * x
    w_nodes = w / math.sqrt(math.pi)
    w_nodes = w_nodes / np.sum(w_nodes)

    p_nodes = np.clip(norm.cdf(y_nodes), 1e-15, 1 - 1e-15)

    h = np.arange(hmax + 1, dtype=int)
    mix = np.zeros(hmax + 1, dtype=float)

    for pj, wj in zip(p_nodes, w_nodes):
        if model == "LD":
            lp = np.array(
                [log_pmf_LD_1pt(int(nbar), int(xx), float(pj), float(q_param)) for xx in h],
                dtype=float
            )
        else:
            lp = np.array(
                [log_pmf_Torri(int(nbar), int(xx), float(pj), float(u), float(v)) for xx in h],
                dtype=float
            )

        pmf = normalize_logpmf(lp)
        mix += wj * pmf

    mix = np.maximum(mix, 0.0)
    mix /= mix.sum()
    return mix




def model_pmf_fixed_n_hierarchical_full(model, key, nbar):
    return model_pmf_fixed_n_hierarchical(model, key, nbar, nbar)




def subset_df_by_period(df, start_year, end_year):
    out = df[(df["Year"] >= start_year) & (df["Year"] <= end_year)].copy()
    out = out.reset_index(drop=True)
    return out




def compute_nbar_for_period(df_sub, classes=("ALL", "SG", "IG")):
    nbar = {}
    for key in classes:
        vals = pd.to_numeric(df_sub[key], errors="coerce").dropna().to_numpy()
        nbar[key] = int(np.rint(np.mean(vals)))
    return nbar




def build_iid_fit_dicts(period_fit_results):
    fits_torri = {}
    fits_ld = {}
    fits_vas = {}
    for key in CLASSES:
        fits_torri[key] = {
            "p": float(period_fit_results["iid"]["Torri"][key]["p"]),
            "u": float(period_fit_results["iid"]["Torri"][key]["u"]),
            "v": float(period_fit_results["iid"]["Torri"][key]["v"]),
        }
        fits_ld[key] = {
            "p": float(period_fit_results["iid"]["LD"][key]["p"]),
            "q": float(period_fit_results["iid"]["LD"][key]["q"]),
        }
        fits_vas[key] = {
            "p": float(period_fit_results["iid"]["Vasicek"][key]["p"]),
            "rho_A": float(period_fit_results["iid"]["Vasicek"][key]["rho_A"]),
        }
    return fits_torri, fits_ld, fits_vas




def build_hier_fit_dicts(period_fit_results):
    fits_hier_ld = {}
    fits_hier_torri = {}
    for key in CLASSES:
        fits_hier_ld[key] = {
            "mu": float(period_fit_results["hier"]["LD"][key]["mu"]),
            "sigma": float(period_fit_results["hier"]["LD"][key]["sigma"]),
            "q": float(period_fit_results["hier"]["LD"][key]["q"]),
            "GH_N": int(period_fit_results["hier"]["LD"][key]["GH_N"]),
        }
        fits_hier_torri[key] = {
            "mu": float(period_fit_results["hier"]["Torri"][key]["mu"]),
            "sigma": float(period_fit_results["hier"]["Torri"][key]["sigma"]),
            "u": float(period_fit_results["hier"]["Torri"][key]["u"]),
            "v": float(period_fit_results["hier"]["Torri"][key]["v"]),
            "GH_N": int(period_fit_results["hier"]["Torri"][key]["GH_N"]),
        }
    return fits_hier_ld, fits_hier_torri




def run_iid_fits_for_period(df_sub):
    out = {"Torri": {}, "LD": {}, "Vasicek": {}}

    for key in CLASSES:
        out["Torri"][key] = fit_one_Torri(df_sub[key], df_sub["D_" + key], maxiter=MAXITER_IID)
        out["LD"][key] = fit_one_LD(df_sub[key], df_sub["D_" + key])
        out["Vasicek"][key] = fit_one_Vas(df_sub[key], df_sub["D_" + key], maxiter=MAXITER_IID)

    rows = []
    for model in ["Torri", "LD", "Vasicek"]:
        for key in CLASSES:
            fit = out[model][key]
            row = {
                "model": model,
                "class": key,
                "nll": float(fit["nll"]),
                "success": bool(fit["success"]),
                "nit": int(fit["nit"]),
                "message": str(fit["message"]),
            }
            for par in ["p", "u", "v", "q", "rho_A", "GH_N"]:
                if par in fit:
                    row[par] = fit[par]
            rows.append(row)

    df_iid = pd.DataFrame(rows).sort_values(["class", "nll"]).reset_index(drop=True)
    return out, df_iid




def run_hier_fits_for_period(df_sub, GH_N=GH_N_HIER):
    out = {"LD": {}, "Torri": {}}

    for key in CLASSES:
        out["LD"][key] = fit_one_hier_LD(
            df_sub[key], df_sub["D_" + key],
            GH_N=GH_N, maxiter=MAXITER_HIER
        )
        out["Torri"][key] = fit_one_hier_Torri(
            df_sub[key], df_sub["D_" + key],
            GH_N=GH_N, maxiter=MAXITER_HIER
        )

    rows = []
    for model in ["LD", "Torri"]:
        for key in CLASSES:
            fit = out[model][key]
            row = {
                "model": f"Hier-{model}",
                "class": key,
                "nll": float(fit["nll"]),
                "success": bool(fit["success"]),
                "nit": int(fit["nit"]),
                "message": str(fit["message"]),
            }
            for par in ["mu", "sigma", "q", "u", "v", "GH_N"]:
                if par in fit:
                    row[par] = fit[par]
            rows.append(row)

    df_hier = pd.DataFrame(rows).sort_values(["class", "nll"]).reset_index(drop=True)
    return out, df_hier




def variance_decomp_vas_fixed_nbar(df, NBAR, fits_vas, keys=("ALL", "SG", "IG")):
    rows = []

    for key in keys:
        nbar = int(NBAR[key])

        n_vals, L_vals = _clean_nL(df, key)
        L_scaled = scaled_counts(n_vals, L_vals, nbar)

        var_data_scaled = float(np.var(L_scaled, ddof=1)) if len(L_scaled) >= 2 else np.nan
        m_hat_data = float(np.sum(L_vals) / np.sum(n_vals))
        var_iid_data_scaled = float(nbar * m_hat_data * (1.0 - m_hat_data))
        var_corr_data_scaled = (
            var_data_scaled - var_iid_data_scaled
            if np.isfinite(var_data_scaled) else np.nan
        )

        p = float(fits_vas[key]["p"])
        rho_A = float(fits_vas[key]["rho_A"])

        a = norm.ppf(np.clip(p, 1e-15, 1.0 - 1e-15))
        q_bar = float(phi2_aa(a, rho_A))
        rho_bar = aggregate_m_q_to_rho(p, q_bar)

        var_iid_model = float(nbar * p * (1.0 - p))
        var_factor_model = float(nbar * (nbar - 1) * (q_bar - p**2))
        var_model_total = var_iid_model + var_factor_model

        rows.append({
            "class": key,
            "nbar": nbar,
            "p": p,
            "rho_A": rho_A,
            "Var_data_scaled": var_data_scaled,
            "m_hat_data": m_hat_data,
            "Var_iid_data_scaled": var_iid_data_scaled,
            "Var_corr_data_scaled": var_corr_data_scaled,
            "m_bar_model": p,
            "q_bar_model": q_bar,
            "rho_bar_model": rho_bar,
            "Var_iid_model": var_iid_model,
            "Var_factor_model": var_factor_model,
            "Var_model_total": var_model_total,
            "Share_iid_within_model": (var_iid_model / var_model_total) if var_model_total > 0 else np.nan,
            "Share_factor_within_model": (var_factor_model / var_model_total) if var_model_total > 0 else np.nan,
            "Iid_over_VarData": (var_iid_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
            "Factor_over_VarData": (var_factor_model / var_data_scaled) if (np.isfinite(var_data_scaled) and var_data_scaled > 0) else np.nan,
        })

    return pd.DataFrame(rows)




def build_compare_decomp_table_with_vas(out_var_hier_ld, out_var_hier_torri, out_var_vas):
    cols_common = ["class", "Iid_over_VarData"]

    df_ld = out_var_hier_ld[[
        "class",
        "Share_iid_within_model",
        "Share_infect_within_model",
        "Share_pt_within_model",
        "Iid_over_VarData",
        "Infect_over_VarData",
        "Pt_over_VarData",
    ]].copy()
    df_ld["model"] = "Hier-Lo–Davis"

    df_t = out_var_hier_torri[[
        "class",
        "Share_iid_within_model",
        "Share_infect_within_model",
        "Share_pt_within_model",
        "Iid_over_VarData",
        "Infect_over_VarData",
        "Pt_over_VarData",
    ]].copy()
    df_t["model"] = "Hier-Torri"

    df_v = out_var_vas[[
        "class",
        "Share_iid_within_model",
        "Share_factor_within_model",
        "Iid_over_VarData",
        "Factor_over_VarData",
    ]].copy()
    df_v["model"] = "Vasicek"

    return pd.concat([df_ld, df_t, df_v], axis=0, ignore_index=True).sort_values(["class", "model"]).reset_index(drop=True)




def pmf_iid_ld_fixed_n(nbar, p, q_param):
    h = np.arange(nbar + 1, dtype=int)
    lp = np.array(
        [log_pmf_LD_1pt(int(nbar), int(x), float(p), float(q_param)) for x in h],
        dtype=float
    )
    return normalize_logpmf(lp)




def pmf_iid_torri_fixed_n(nbar, p, u, v):
    h = np.arange(nbar + 1, dtype=int)
    lp = np.array(
        [log_pmf_Torri(int(nbar), int(x), float(p), float(u), float(v)) for x in h],
        dtype=float
    )
    return normalize_logpmf(lp)




def compute_tail_table_for_period(
    df_sub,
    iid_fit_results,
    hier_fit_results,
    NBAR,
    alphas=(0.95, 0.99),
):
    fits_torri, fits_ld, fits_vas = build_iid_fit_dicts({"iid": iid_fit_results})
    fits_hier_ld, fits_hier_torri = build_hier_fit_dicts({"hier": hier_fit_results})

    rows = []

    for key in CLASSES:
        nbar = int(NBAR[key])

        n_vals = df_sub[key].astype(int).values
        L_vals = df_sub["D_" + key].astype(int).values
        L_scaled = scaled_counts(n_vals, L_vals, nbar)

        _, pmf_data = empirical_pmf_from_counts(L_scaled)
        pmf_data_full = np.zeros(nbar + 1, dtype=float)
        pmf_data_full[:len(pmf_data)] = pmf_data

        pmf_v = pmf_vas_fixed_n(nbar, fits_vas[key]["p"], fits_vas[key]["rho_A"])
        pmf_iid_ld = pmf_iid_ld_fixed_n(nbar, fits_ld[key]["p"], fits_ld[key]["q"])
        pmf_iid_t = pmf_iid_torri_fixed_n(nbar, fits_torri[key]["p"], fits_torri[key]["u"], fits_torri[key]["v"])
        pmf_hld = pmf_hier_ld_fixed_n(
            nbar,
            fits_hier_ld[key]["mu"],
            fits_hier_ld[key]["sigma"],
            fits_hier_ld[key]["q"],
            fits_hier_ld[key]["GH_N"],
        )
        pmf_ht = pmf_hier_torri_fixed_n(
            nbar,
            fits_hier_torri[key]["mu"],
            fits_hier_torri[key]["sigma"],
            fits_hier_torri[key]["u"],
            fits_hier_torri[key]["v"],
            fits_hier_torri[key]["GH_N"],
        )

        model_pmfs = [
            ("Data(scaled)", pmf_data_full),
            ("IID-Vasicek", pmf_v),
            ("IID-Lo–Davis", pmf_iid_ld),
            ("IID-Torri", pmf_iid_t),
            ("Hier-Lo–Davis", pmf_hld),
            ("Hier-Torri", pmf_ht),
        ]

        for alpha in alphas:
            for model_name, pmf in model_pmfs:
                met = tail_metrics_from_pmf(pmf, alpha=alpha, strict_es=False)
                rows.append({
                    "class": key,
                    "alpha": alpha,
                    "model": model_name,
                    "VaR": met["VaR"],
                    "ES": met["ES"],
                    "ES_minus_VaR": met["ES_minus_VaR"],
                    "ES_over_VaR": met["ES_over_VaR"],
                })

    return pd.DataFrame(rows).sort_values(["class", "alpha", "model"]).reset_index(drop=True)




def run_period_analysis(df, period_name, start_year, end_year):
    df_sub = subset_df_by_period(df, start_year, end_year)
    NBAR = compute_nbar_for_period(df_sub, CLASSES)

    print(f"\n===== {period_name}: {start_year}-{end_year} =====")
    print(f"N_years = {len(df_sub)}")
    print(f"NBAR = {NBAR}")

    # IID
    iid_fit_results, df_iid = run_iid_fits_for_period(df_sub)

    # Hierarchical
    hier_fit_results, df_hier = run_hier_fits_for_period(df_sub, GH_N=GH_N_HIER)

    # Build dicts
    fits_torri, fits_ld, fits_vas = build_iid_fit_dicts({"iid": iid_fit_results})
    fits_hier_ld, fits_hier_torri = build_hier_fit_dicts({"hier": hier_fit_results})

    # Variance decomposition
    out_var_hier_ld = variance_decomp_hier_ld_fixed_nbar(
        df_sub, NBAR, fits_hier_ld, CLASSES
    )

    out_var_hier_torri = variance_decomp_hier_torri_fixed_nbar(
        df_sub, NBAR, fits_hier_torri, CLASSES
    )

    out_var_vas = variance_decomp_vas_fixed_nbar(
        df_sub, NBAR, fits_vas, CLASSES
    )

    df_var_compare = build_compare_decomp_table_with_vas(
        out_var_hier_ld,
        out_var_hier_torri,
        out_var_vas,
    )

    # Tail table
    df_tail = compute_tail_table_for_period(
        df_sub=df_sub,
        iid_fit_results=iid_fit_results,
        hier_fit_results=hier_fit_results,
        NBAR=NBAR,
        alphas=ALPHAS,
    )

    return {
        "df_sub": df_sub,
        "NBAR": NBAR,
        "iid_fit_results": iid_fit_results,
        "hier_fit_results": hier_fit_results,
        "iid_fit_table": df_iid,
        "hier_fit_table": df_hier,
        "var_hier_ld": out_var_hier_ld,
        "var_hier_torri": out_var_hier_torri,
        "var_vas": out_var_vas,
        "var_compare": df_var_compare,
        "tail_table": df_tail,
    }




def compact_nll_summary(period_results):
    rows = []
    for pname, res in period_results.items():
        tmp1 = res["iid_fit_table"][["class", "model", "nll"]].copy()
        tmp1["period"] = pname
        tmp1["type"] = "IID"

        tmp2 = res["hier_fit_table"][["class", "model", "nll"]].copy()
        tmp2["period"] = pname
        tmp2["type"] = "Hier"

        rows.append(tmp1)
        rows.append(tmp2)

    out = pd.concat(rows, axis=0, ignore_index=True)
    return out[["period", "type", "class", "model", "nll"]].sort_values(
        ["period", "type", "class", "nll"]
    ).reset_index(drop=True)




def prepare_plot_data(df_cmp_decomp):
    classes = ["ALL", "SG", "IG"]
    models = ["Lo–Davis", "Torri"]

    labels = []
    iid_vals = []
    inf_vals = []
    pt_vals = []

    for c in classes:
        for m in models:
            sub = df_cmp_decomp[
                (df_cmp_decomp["class"] == c) &
                (df_cmp_decomp["model"] == m)
            ]

            labels.append(f"{c}\n{'LD' if m == 'Lo–Davis' else 'Torri'}")

            iid_vals.append(sub["Iid_over_VarData"].values[0])
            inf_vals.append(sub["Infect_over_VarData"].values[0])
            pt_vals.append(sub["Pt_over_VarData"].values[0])

    return (
        labels,
        np.array(iid_vals),
        np.array(inf_vals),
        np.array(pt_vals)
    )




def to_float_dict(d):
    out = {}
    for k, v in d.items():
        try:
            out[k] = float(v)
        except Exception:
            out[k] = v
    return out


