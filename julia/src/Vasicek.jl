module Vasicek

using StatsFuns: normcdf, norminvcdf
using LogExpFunctions: logsumexp
using FastGaussQuadrature: gausshermite
using ..ProbitNormal: log_comb

export vasicek_to_probit_normal, probit_normal_to_vasicek, log_pmf_Vas_1pt, nll_vasicek

"""
    vasicek_to_probit_normal(p0::Float64, rho_A::Float64)::Tuple{Float64, Float64}

Convert Vasicek parameters (p0, rho_A) to Probit-Normal (mu, sigma) parameters.
- mu = norminvcdf(p0) / sqrt(1 - rho_A)
- sigma = sqrt(rho_A / (1 - rho_A))
"""
function vasicek_to_probit_normal(p0::Float64, rho_A::Float64)::Tuple{Float64, Float64}
    if !(0.0 < p0 < 1.0) || !(0.0 <= rho_A < 1.0)
        throw(ArgumentError("Require 0 < p0 < 1 and 0 <= rho_A < 1"))
    end
    val_inv = sqrt(1.0 - rho_A)
    mu = norminvcdf(p0) / val_inv
    sigma = sqrt(rho_A / (1.0 - rho_A))
    return (mu, sigma)
end

"""
    probit_normal_to_vasicek(mu::Float64, sigma::Float64)::Tuple{Float64, Float64}

Convert Probit-Normal parameters (mu, sigma) to Vasicek (p0, rho_A) parameters.
- p0 = normcdf(mu / sqrt(1 + sigma^2))
- rho_A = sigma^2 / (1 + sigma^2)
"""
function probit_normal_to_vasicek(mu::Float64, sigma::Float64)::Tuple{Float64, Float64}
    if sigma < 0.0
        throw(ArgumentError("sigma must be >= 0"))
    end
    var_y = sigma^2
    p0 = normcdf(mu / sqrt(1.0 + var_y))
    rho_A = var_y / (1.0 + var_y)
    return (p0, rho_A)
end

"""
    log_pmf_Vas_1pt(n::Int, h::Int, p::Float64, rho_A::Float64; GH_N::Int=100)::Float64

Exact single-observation Vasicek log PMF matching Python `models.log_pmf_Vas_1pt`.
"""
function log_pmf_Vas_1pt(n::Int, h::Int, p::Float64, rho_A::Float64; GH_N::Int=100)::Float64
    if h < 0 || h > n
        return -Inf
    end

    eps_val = 1e-15
    p_clamped = clamp(p, eps_val, 1.0 - eps_val)
    rho_clamped = clamp(rho_A, eps_val, 1.0 - 1e-12)

    a = norminvcdf(p_clamped)
    sr = sqrt(rho_clamped)
    tr = sqrt(1.0 - rho_clamped)

    gh_x, gh_w = gausshermite(GH_N)
    
    logC = log_comb(n, h)
    K = length(gh_x)
    log_terms = Vector{Float64}(undef, K)

    @inbounds for i in 1:K
        z = sqrt(2.0) * gh_x[i]
        pi_val = normcdf((a - sr * z) / tr)
        pi_clamped = clamp(pi_val, eps_val, 1.0 - eps_val)
        
        log_bin = logC + h * log(pi_clamped) + (n - h) * log1p(-pi_clamped)
        log_terms[i] = log(gh_w[i]) + log_bin
    end

    return logsumexp(log_terms) - 0.5 * log(pi)
end

"""
    nll_vasicek(p0::Float64, rho_A::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64

Compute Negative Log-Likelihood for Vasicek model over grouped annual data.
"""
function nll_vasicek(p0::Float64, rho_A::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64
    tot_nll = 0.0
    for i in 1:length(n_arr)
        lp = log_pmf_Vas_1pt(n_arr[i], h_arr[i], p0, rho_A; GH_N=GH_N)
        if !isfinite(lp)
            return 1e100
        end
        tot_nll -= c_arr[i] * lp
    end
    return tot_nll
end

end # module Vasicek
