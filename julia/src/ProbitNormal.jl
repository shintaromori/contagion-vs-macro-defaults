module ProbitNormal

using StatsFuns: normcdf, loggamma
using LogExpFunctions: logsumexp
using ..Quadrature: get_hermite_quadrature

export log_comb, nll_probit_normal, nll_probit_normal_unconstrained

"""
    log_comb(n::Int, h::Int)::Float64

Compute log binomial coefficient log(n! / (h! * (n-h)!)).
"""
function log_comb(n::Int, h::Int)::Float64
    if h < 0 || h > n
        return -Inf
    end
    return loggamma(n + 1) - loggamma(h + 1) - loggamma(n - h + 1)
end

"""
    nll_probit_normal(mu::Float64, sigma::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64

Compute Negative Log-Likelihood (NLL) for single Probit-Normal NULL model (M0) with physical parameters (mu, sigma).
"""
function nll_probit_normal(mu::Float64, sigma::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64
    if sigma <= 0.0
        return Inf
    end

    x, log_w = get_hermite_quadrature(GH_N)
    
    # y_nodes = mu + sqrt(2) * sigma * x
    sqrt2_sigma = sqrt(2.0) * sigma
    y_nodes = mu .+ sqrt2_sigma .* x
    
    # Clip CDF values to [1e-15, 1 - 1e-15] matching Python
    eps_val = 1e-15
    p_nodes = clamp.(normcdf.(y_nodes), eps_val, 1.0 - eps_val)
    
    log_p = log.(p_nodes)
    log_1mp = log1p.(-p_nodes)
    
    tot_nll = 0.0
    N = length(n_arr)
    
    # Buffer for node terms
    K = length(x)
    node_terms = Vector{Float64}(undef, K)
    
    for i in 1:N
        n_i = n_arr[i]
        h_i = h_arr[i]
        cnt_i = c_arr[i]
        
        log_c = log_comb(n_i, h_i)
        
        # log_terms = log_c + h * log_p + (n - h) * log_1mp
        @inbounds for k in 1:K
            node_terms[k] = log_w[k] + log_c + h_i * log_p[k] + (n_i - h_i) * log_1mp[k]
        end
        
        log_prob = logsumexp(node_terms)
        tot_nll -= cnt_i * log_prob
    end
    
    return tot_nll
end

"""
    nll_probit_normal_unconstrained(params::Vector{Float64}, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64

Unconstrained parameterization wrapper: params = [mu, log_sigma].
"""
function nll_probit_normal_unconstrained(params::Vector{Float64}, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64
    mu = params[1]
    sigma = exp(params[2])
    return nll_probit_normal(mu, sigma, n_arr, h_arr, c_arr; GH_N=GH_N)
end

end # module ProbitNormal
