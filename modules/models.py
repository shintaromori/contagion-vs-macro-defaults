"""Core probability models and fitting routines for the IDM project."""

import math
import numpy as np

import pandas as pd
from math import lgamma
from scipy.optimize import minimize, brentq
from scipy.stats import norm, multivariate_normal
from numpy.polynomial.hermite import hermgauss


def log_comb(n: int, k: int) -> float:
    """log( n choose k )"""
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)




def logsumexp_arr(logw: np.ndarray) -> float:
    """stable log(sum(exp(logw))) for 1D array"""
    m = np.max(logw)
    if np.isneginf(m):
        return -np.inf
    return m + np.log(np.sum(np.exp(logw - m)))




def logsumexp2(x: float, y: float) -> float:
    """log(exp(x)+exp(y))"""
    m = max(x, y)
    if np.isneginf(m):
        return -np.inf
    return m + np.log(np.exp(x-m) + np.exp(y-m))




def logsubexp(a: float, b: float) -> float:
    """log(exp(a)-exp(b)) with a>b; otherwise -inf"""
    if b >= a:
        return -np.inf
    return a + np.log1p(-np.exp(b - a))




def normalize_logpmf(logpmf: np.ndarray) -> np.ndarray:
    """Normalize log-PMF array to PMF (handles -inf safely)."""
    logpmf = np.where(np.isfinite(logpmf), logpmf, -np.inf)
    m = np.max(logpmf)
    if not np.isfinite(m):
        return np.full_like(logpmf, np.nan, dtype=float)
    z = np.sum(np.exp(logpmf - m))
    return np.exp(logpmf - (m + np.log(z)))




def var_from_pmf(pmf: np.ndarray, alpha: float = 0.99, side: str = "left") -> int:
    """VaRα = inf{h: P(L<=h) >= α}."""
    cdf = np.cumsum(pmf)
    return int(np.searchsorted(cdf, alpha, side=side))




def es_from_pmf(pmf: np.ndarray, var_h: int, strict: bool = False) -> float:
    """
    ESα: default is E[L | L >= VaRα] (strict=False).
    strict=True -> E[L | L > VaRα]
    """
    pmf = np.asarray(pmf, dtype=float)
    x = np.arange(len(pmf))
    start = var_h + 1 if strict else var_h
    tail = pmf[start:]
    tp = tail.sum()
    if tp <= 0:
        return np.nan
    return float((x[start:] * tail).sum() / tp)




def moments_from_pmf(pmf: np.ndarray):
    x = np.arange(len(pmf))
    mean = float((x*pmf).sum())
    var = float(((x-mean)**2 * pmf).sum())
    return mean, var




def m_rho_from_pmf(pmf: np.ndarray):
    """Return (m, rho, Var(L)) where m=E[L]/n, rho is pairwise default corr (exchangeable proxy)."""
    pmf = np.asarray(pmf, dtype=float)
    s = pmf.sum()
    if s <= 0:
        return np.nan, np.nan, np.nan
    pmf = pmf / s

    n = len(pmf) - 1
    x = np.arange(n+1)

    EL  = float((x * pmf).sum())
    EL2 = float(((x**2) * pmf).sum())
    varL = EL2 - EL**2

    if n <= 0:
        return np.nan, np.nan, varL

    m = EL / n
    denom = n*(n-1)*m*(1-m)
    if denom <= 1e-15:
        return m, np.nan, varL

    rho = (varL - n*m*(1-m)) / denom
    return m, rho, varL




def survival_from_pmf(pmf: np.ndarray) -> np.ndarray:
    """S[h] = P(L >= h)"""
    return np.flip(np.cumsum(np.flip(pmf)))




def tail_metrics_from_pmf(pmf: np.ndarray, alpha: float = 0.99, strict_es: bool = False):
    v = var_from_pmf(pmf, alpha=alpha)
    e = es_from_pmf(pmf, v, strict=strict_es)
    return dict(
        VaR=v,
        ES=e,
        ES_minus_VaR=(e - v) if np.isfinite(e) else np.nan,
        ES_over_VaR=(e / v) if (np.isfinite(e) and v > 0) else np.nan
    )




def pmf_Torri(n: int, p: float, u: float, v: float) -> np.ndarray:
    """
    Torri (infectious default) PMF for L_n.
    Numerically stable (log-sum-exp).
    """
    a = p*(1-v)
    b = (1-p)*(1-u)
    c = (1-p)*u
    bc = b + c  # = 1-p

    pmf_log = np.full(n+1, -np.inf, dtype=float)

    log_bc = -np.inf if bc <= 0 else np.log(bc)

    for h in range(n+1):
        # term0 = a^h * (b+c)^(n-h)
        if h == 0:
            log_a_h = 0.0
        else:
            log_a_h = -np.inf if a <= 0 else h*np.log(a)
        log_term0 = log_a_h + (n-h)*log_bc

        # term1 = ((p+b)^h - (a+b)^h) * c^(n-h)
        if h == 0:
            log_term1 = -np.inf
        else:
            # c^(n-h)
            if (n-h) == 0:
                log_c_pow = 0.0
            else:
                log_c_pow = -np.inf if c <= 0 else (n-h)*np.log(c)

            X = p + b
            Y = a + b

            if log_c_pow == -np.inf or X <= 0 or Y <= 0:
                log_term1 = -np.inf
            else:
                logXh = h*np.log(X)
                logYh = h*np.log(Y)
                log_diff = logsubexp(logXh, logYh)
                log_term1 = -np.inf if log_diff == -np.inf else (log_diff + log_c_pow)

        pmf_log[h] = log_comb(n, h) + logsumexp2(log_term0, log_term1)

    return normalize_logpmf(pmf_log)




def moments_from_params_Torri(n: int, p: float, u: float, v: float):
    """Closed-form (m, rho) for Torri parameterization used in this notebook."""
    b = (1-p)*(1-u)
    s = 1 - p*v

    # mean (fraction)
    m = p + b*(1 - s**(n-1))

    # pair moment term
    a_alt = s + p - 1
    s_pow = np.exp((n-2)*np.log(s)) if n > 2 else 1.0
    E = (p+b)**2 - b*(2*a_alt + b)*s_pow

    cov = E - m**2
    den = m*(1-m)
    rho = np.nan if den <= 1e-15 else cov/den
    return float(m), float(rho)




def find_params_on_iso_Torri(n: int, m_target: float, rho_target: float, p_grid: np.ndarray):
    """
    Solve (u,v) on iso-(m,rho) manifold by scanning p and solving for s.
    Returns list of dicts with keys: p,u,v,s,b
    """
    sols = []

    def bisect_root(f, lo, hi, it=80):
        flo, fhi = f(lo), f(hi)
        if not (np.isfinite(flo) and np.isfinite(fhi)) or flo*fhi > 0:
            return None
        for _ in range(it):
            mid = 0.5*(lo+hi)
            fmid = f(mid)
            if not np.isfinite(fmid):
                mid = np.nextafter(mid, lo)
                fmid = f(mid)
                if not np.isfinite(fmid):
                    return None
            if flo*fmid <= 0:
                hi, fhi = mid, fmid
            else:
                lo, flo = mid, fmid
        return 0.5*(lo+hi)

    for p in p_grid:
        if not (0 < p < m_target):
            continue

        s_lo = (1 - p) + 1e-10
        s_hi = 1 - 1e-10
        if s_lo >= s_hi:
            continue

        def b_of_s(s):
            denom = 1 - s**(n-1)
            if denom <= 0:
                return np.nan
            return (m_target - p) / denom

        def F(s):
            b = b_of_s(s)
            if not np.isfinite(b):
                return np.nan
            if not (0 <= b <= 1-p + 1e-12):
                return np.nan

            a_alt = s + p - 1
            s_pow = np.exp((n-2)*np.log(s)) if n > 2 else 1.0

            E = (p+b)**2 - b*(2*a_alt + b)*s_pow
            cov = E - (m_target**2)
            den = m_target*(1-m_target)
            rho = cov/den if den > 1e-15 else np.nan
            return rho - rho_target

        S = np.linspace(s_lo, s_hi, 600)
        Fv = np.array([F(s) for s in S])

        for k in range(len(S)-1):
            f1, f2 = Fv[k], Fv[k+1]
            if not (np.isfinite(f1) and np.isfinite(f2)):
                continue
            if f1 == 0:
                s_star = S[k]
            elif f1*f2 < 0:
                s_star = bisect_root(F, S[k], S[k+1])
                if s_star is None:
                    continue
            else:
                continue

            b = b_of_s(s_star)
            u = 1 - b/(1-p)
            v = (1 - s_star)/p
            if not (0 <= u <= 1 and 0 <= v <= 1):
                continue

            sols.append({"p": float(p), "u": float(u), "v": float(v),
                         "s": float(s_star), "b": float(b)})

    sols.sort(key=lambda d: d["p"])
    uniq = []
    for d in sols:
        if len(uniq) == 0:
            uniq.append(d); continue
        e = uniq[-1]
        if abs(d["u"]-e["u"]) + abs(d["v"]-e["v"]) > 5e-4:
            uniq.append(d)
    return uniq




def pmf_LD(n: int, p: float, q: float) -> np.ndarray:
    """
    Lo–Davis Infectious Default Model PMF via mixture:
      K ~ Bin(n,p), L|K=k = k + Bin(n-k, r_k), r_k = 1 - (1-q)^k
    Stable implementation using expm1/log1p and log-sum-exp.
    """
    if not (0.0 <= p <= 1.0 and 0.0 <= q <= 1.0):
        raise ValueError("p and q must be in [0,1].")

    # log P(K=k)
    logP_K = np.full(n+1, -np.inf, dtype=float)
    if p == 0.0:
        logP_K[0] = 0.0
    elif p == 1.0:
        logP_K[n] = 0.0
    else:
        log_p = np.log(p)
        log_1mp = np.log1p(-p)
        for k in range(n+1):
            logP_K[k] = log_comb(n, k) + k*log_p + (n-k)*log_1mp

    # r_k = 1 - (1-q)^k  (stable)
    r = np.zeros(n+1, dtype=float)
    if q == 0.0:
        r[:] = 0.0
    elif q == 1.0:
        r[0] = 0.0
        r[1:] = 1.0
    else:
        l = np.log1p(-q)
        ks = np.arange(n+1, dtype=float)
        r = -np.expm1(ks * l)

    pmf_log = np.full(n+1, -np.inf, dtype=float)

    for h in range(n+1):
        # sum over k=0..h
        terms = np.full(h+1, -np.inf, dtype=float)
        for k in range(h+1):
            if not np.isfinite(logP_K[k]):
                continue
            t = h - k
            nk = n - k
            rk = r[k]
            if t < 0 or t > nk:
                continue

            if rk == 0.0:
                if t != 0:
                    continue
                log_bin = 0.0
            elif rk == 1.0:
                if t != nk:
                    continue
                log_bin = 0.0
            else:
                log_rk = np.log(rk)
                log_1mrk = np.log1p(-rk)
                log_bin = log_comb(nk, t) + t*log_rk + (nk-t)*log_1mrk

            terms[k] = logP_K[k] + log_bin

        pmf_log[h] = logsumexp_arr(terms)

    return normalize_logpmf(pmf_log)




def m_rho_LD(n: int, p: float, q: float):
    """
    Closed-form (m, rhoZ) for Lo–Davis.
    m = 1 - (1-p)*(1-qp)^(n-1)
    Cov = (1-p)^2 * [ (1 - 2pq + p q^2)^(n-2) - (1-pq)^(2(n-1)) ]
    rho = Cov / (m(1-m))
    """
    m = 1.0 - (1.0 - p) * (1.0 - q*p)**(n-1)
    if n < 2:
        return float(m), np.nan
    termA = (1.0 - 2.0*p*q + p*(q**2))**(n-2)
    termB = (1.0 - p*q)**(2*(n-1))
    cov = (1.0 - p)**2 * (termA - termB)
    var = m*(1.0-m)
    rho = cov/var if var > 1e-15 else np.nan
    return float(m), float(rho)




def solve_pq_for_mrho_LD(n: int, m_target: float, rho_target: float,
                         p_bounds=(1e-10, 1-1e-10), q_bounds=(1e-10, 1-1e-10),
                         n_tries: int = 60, w_rho: float = 1.0):
    """Fit (p,q) to match (m_target, rho_target) by multi-start L-BFGS-B."""
    p_lo, p_hi = p_bounds
    q_lo, q_hi = q_bounds

    def objective(x):
        p, q = x
        if not (p_lo <= p <= p_hi and q_lo <= q <= q_hi):
            return 1e9
        m, rho = m_rho_LD(n, p, q)
        if not np.isfinite(m) or not np.isfinite(rho):
            return 1e9
        return (m - m_target)**2 + w_rho*(rho - rho_target)**2

    best = None
    rng = np.random.default_rng(SEED)

    init_list = [
        (min(max(m_target, p_lo), p_hi), 0.05),
        (0.01, 0.10), (0.05, 0.05), (0.10, 0.10), (0.30, 0.05), (0.50, 0.05),
    ]
    for _ in range(max(0, n_tries - len(init_list))):
        init_list.append((rng.uniform(p_lo, p_hi), rng.uniform(q_lo, q_hi)))

    for (p0, q0) in init_list:
        res = minimize(
            objective,
            x0=np.array([p0, q0], dtype=float),
            method="L-BFGS-B",
            bounds=[(p_lo, p_hi), (q_lo, q_hi)],
        )
        if not res.success:
            continue

        p_hat, q_hat = res.x
        m_hat, rho_hat = m_rho_LD(n, p_hat, q_hat)
        err_m = m_hat - m_target
        err_r = rho_hat - rho_target
        score = err_m**2 + w_rho*(err_r**2)

        cand = dict(
            p=float(p_hat), q=float(q_hat),
            m=float(m_hat), rho=float(rho_hat),
            err_m=float(err_m), err_rho=float(err_r),
            score=float(score), fun=float(res.fun),
            success=bool(res.success), message=str(res.message)
        )
        if (best is None) or (cand["score"] < best["score"]):
            best = cand

    if best is None:
        raise RuntimeError("最適化が収束しませんでした。(m, rho) が Lo–Davis で可到達か確認してください。")
    return best




def check_mrho_formula_vs_pmf_LD(n: int, p: float, q: float):
    pmf = pmf_LD(n, p, q)
    m_pmf, rho_pmf, _ = m_rho_from_pmf(pmf)
    m_f, rho_f = m_rho_LD(n, p, q)
    return dict(m_pmf=m_pmf, m_formula=m_f, diff_m=m_pmf-m_f,
                rho_pmf=rho_pmf, rho_formula=rho_f, diff_rho=rho_pmf-rho_f,
                pmf_sum=float(pmf.sum()))




def phi2_aa(a: float, rhoA: float) -> float:
    cov = np.array([[1.0, rhoA],
                    [rhoA, 1.0]])
    return float(multivariate_normal.cdf([a, a], mean=[0.0, 0.0], cov=cov))




def default_corr_from_asset_corr(m: float, rhoA: float) -> float:
    a = norm.ppf(m)
    q = phi2_aa(a, rhoA)
    return float((q - m*m) / (m*(1.0 - m)))




def invert_asset_corr_from_default_corr(m: float, rho_target: float, tol: float = 1e-10) -> float:
    """Solve rhoA in [0, 0.999999] such that Corr(Zi,Zj)=rho_target."""
    if not (0.0 < m < 1.0):
        raise ValueError("m must be in (0,1).")
    if not (-1.0 < rho_target < 1.0):
        raise ValueError("rho_target must be in (-1,1).")

    def F(rhoA):
        return default_corr_from_asset_corr(m, rhoA) - rho_target

    return float(brentq(F, 0.0, 0.999999, xtol=tol, rtol=tol, maxiter=200))




def pmf_Vas(n: int, m: float, rhoA: float, n_gh: int = 60) -> np.ndarray:
    """
    L|Y=y ~ Bin(n, p(y)), p(y)=Phi((a - sqrt(rhoA)*y)/sqrt(1-rhoA))
    PMF = E[ BinPMF(h; n, p(Y)) ] using Gauss-Hermite + log-sum-exp.
    """
    a = norm.ppf(m)
    s = np.sqrt(rhoA)
    t = np.sqrt(1.0 - rhoA)

    x, w = hermgauss(n_gh)
    y = np.sqrt(2.0) * x
    log_w = np.log(w) - 0.5*np.log(np.pi)

    p = norm.cdf((a - s*y) / t)
    p = np.clip(p, 1e-15, 1.0 - 1e-15)
    logp = np.log(p)
    log1mp = np.log1p(-p)

    hs = np.arange(n+1)
    logC = np.array([log_comb(n, h) for h in hs], dtype=float)

    log_bin = logC[None, :] + hs[None, :]*logp[:, None] + (n-hs)[None, :]*log1mp[:, None]
    log_terms = log_w[:, None] + log_bin

    m0 = np.max(log_terms, axis=0)
    pmf = np.exp(m0 + np.log(np.sum(np.exp(log_terms - m0[None, :]), axis=0)))

    pmf = np.maximum(pmf, 0.0)
    pmf = pmf / pmf.sum()
    return pmf




def log_pmf_Torri(n: int, h: int, p: float, u: float, v: float) -> float:
    # constraints (p strictly inside to avoid log(0) issues in bc=1-p)
    if not (0.0 < p < 1.0 and 0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        return -np.inf
    if not (0 <= h <= n):
        return -np.inf

    a = p*(1.0 - v)
    b = (1.0 - p)*(1.0 - u)
    c = (1.0 - p)*u
    bc = b + c  # = 1-p

    # --- term0 = a^h * (b+c)^(n-h) ---
    # log term0
    if h > 0:
        if a <= 0.0:
            t0 = -np.inf
        else:
            t0 = h*math.log(a) + (n-h)*math.log(bc)
    else:
        # a^0 = 1
        t0 = (n-h)*math.log(bc)

    # --- term1 = ((p+b)^h - (a+b)^h) * c^(n-h) ---
    pb = p + b
    ab = a + b

    if h == 0:
        t1 = -np.inf
    else:
        # if v=0 -> a=p -> pb==ab -> difference 0
        if pb <= 0.0 or ab <= 0.0:
            t1 = -np.inf
        else:
            lp = h*math.log(pb)
            la = h*math.log(ab)

            # stable log(exp(lp)-exp(la)) (requires lp>la)
            diff_log = logsubexp(lp, la)
            if diff_log == -np.inf:
                t1 = -np.inf
            else:
                if (n-h) > 0:
                    if c <= 0.0:
                        t1 = -np.inf
                    else:
                        t1 = diff_log + (n-h)*math.log(c)
                else:
                    t1 = diff_log

    # combine
    inner = logsumexp(t0, t1)
    if inner == -np.inf:
        return -np.inf

    return log_comb(n, h) + inner




def logsumexp_vec(logw):
    m = np.max(logw)
    if np.isneginf(m):
        return -np.inf
    return float(m + np.log(np.sum(np.exp(logw - m))))




def log_pmf_LD_1pt(n, h, p, q):
    """
    LD mixture:
      K ~ Bin(n,p)
      L|K=k = k + Bin(n-k, r_k),  r_k = 1 - (1-q)^k
    Return log P(L=h) by summing over k=0..h
    """
    if not (0 <= h <= n):
        return -np.inf

    # numerical guards
    eps = 1e-15
    p = float(np.clip(p, eps, 1.0 - eps))
    q = float(np.clip(q, eps, 1.0 - eps))

    log_p  = math.log(p)
    log_1p = math.log1p(-p)
    log_1q = math.log1p(-q)

    terms = np.full(h+1, -np.inf, dtype=float)

    for k in range(0, h+1):
        t = h - k
        nk = n - k
        if t < 0 or t > nk:
            continue

        # log P(K=k)
        log_pk = log_comb(n, k) + k*log_p + (n-k)*log_1p

        # r_k = 1 - (1-q)^k  (stable)
        if k == 0:
            if t != 0:
                continue
            log_bin = 0.0
        else:
            rk = -math.expm1(k * log_1q)  # stable 1-(1-q)^k
            # handle extreme rk
            if rk <= 0.0:
                if t != 0:
                    continue
                log_bin = 0.0
            elif rk >= 1.0:
                if t != nk:
                    continue
                log_bin = 0.0
            else:
                log_r  = math.log(rk)
                log_1r = math.log1p(-rk)
                log_bin = log_comb(nk, t) + t*log_r + (nk-t)*log_1r

        terms[k] = log_pk + log_bin

    return logsumexp_vec(terms)





def log_pmf_Vas_1pt(n: int, h: int, p: float, rho_A: float, GH_N: int = 100) -> float:
    if not (0 <= h <= n):
        return -np.inf

    eps = 1e-15
    p = float(np.clip(p, eps, 1.0 - eps))
    rho_A = float(np.clip(rho_A, eps, 1.0 - 1e-12))

    a  = norm.ppf(p)
    sr = math.sqrt(rho_A)
    tr = math.sqrt(1.0 - rho_A)

    gh_x, gh_w = hermgauss(int(GH_N))
    mask = np.isfinite(gh_w) & (gh_w > 0) & np.isfinite(gh_x)
    x = gh_x[mask]
    w = gh_w[mask]
    if len(w) == 0:
        return -np.inf

    logC = log_comb(n, h)
    log_terms = np.empty(len(w), dtype=float)

    for i, (xi, wi) in enumerate(zip(x, w)):
        z  = math.sqrt(2.0) * float(xi)
        pi = norm.cdf((a - sr*z) / tr)
        pi = float(np.clip(pi, eps, 1.0 - eps))
        log_bin = logC + h*math.log(pi) + (n-h)*math.log1p(-pi)
        log_terms[i] = math.log(float(wi)) + log_bin

    return logsumexp_vec(log_terms) - 0.5*math.log(math.pi)



def fit_one_Vas(series_n, series_h, x0=None, maxiter=8000, GH_N: int = 100):
    tmp = pd.DataFrame({
        "n": series_n.to_numpy(dtype=int),
        "h": series_h.to_numpy(dtype=int),
    })
    g = tmp.value_counts().reset_index(name="cnt")
    n_arr = g["n"].to_numpy(dtype=int)
    h_arr = g["h"].to_numpy(dtype=int)
    c_arr = g["cnt"].to_numpy(dtype=int)

    def neg_ll_vas(theta):
        p     = 1.0/(1.0 + np.exp(-theta[0]))
        rho_A = 1.0/(1.0 + np.exp(-theta[1]))

        ll = 0.0
        for n, h, c in zip(n_arr, h_arr, c_arr):
            lp = log_pmf_Vas_1pt(int(n), int(h), float(p), float(rho_A), GH_N=GH_N)
            if not np.isfinite(lp):
                return 1e100
            ll += int(c) * lp
        return -ll

    if x0 is None:
        x0 = np.array([math.log(0.01/0.99), math.log(0.05/0.95)], dtype=float)

    res = minimize(neg_ll_vas, x0, method="Nelder-Mead", options={"maxiter": maxiter})

    th = res.x
    p     = 1.0/(1.0 + np.exp(-th[0]))
    rho_A = 1.0/(1.0 + np.exp(-th[1]))

    return {
        "p": float(p),
        "rho_A": float(rho_A),
        "nll": float(res.fun),
        "success": bool(res.success),
        "nit": int(res.nit),
        "message": str(res.message),
        "GH_N": int(GH_N),
    }
