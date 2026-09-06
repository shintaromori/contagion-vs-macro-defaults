module GaussianMixture

using StatsFuns: normcdf
using LogExpFunctions: logsumexp, logaddexp, softplus, logistic, logit
using ..Quadrature: get_hermite_quadrature
using ..ProbitNormal: log_comb

export nll_probit_normal_mixture_physical, unconstrained_to_mixture_params, nll_probit_normal_mixture_unconstrained
export inv_softplus, mixture_to_unconstrained_params

"""
    inv_softplus(y::Float64)::Float64

Numerically stable inverse of softplus(x) = log(1 + exp(x)).
For y > 20.0, softplus(y) ≈ y.
For small y > 0, log(expm1(y)) provides exact precision without catastrophic cancellation.
"""
function inv_softplus(y::Float64)::Float64
    if y <= 0.0
        throw(ArgumentError("Require y > 0 for inv_softplus"))
    elseif y > 20.0
        return y
    else
        return log(expm1(y))
    end
end

"""
    nll_probit_normal_mixture_physical(mu1::Float64, mu2::Float64, sigma1::Float64, sigma2::Float64, pi_val::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64

Direct physical Negative Log-Likelihood evaluation for 2-component Gaussian mixture environmental NULL model (M3).
y_t ~ pi N(mu1, sigma1^2) + (1-pi) N(mu2, sigma2^2)
p_t = Phi(y_t)
L_t | p_t, n_t ~ Binomial(n_t, p_t)
"""
function nll_probit_normal_mixture_physical(mu1::Float64, mu2::Float64, sigma1::Float64, sigma2::Float64, pi_val::Float64, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64
    if !(0.0 < pi_val < 1.0) || sigma1 <= 0.0 || sigma2 <= 0.0
        return Inf
    end

    x, log_w = get_hermite_quadrature(GH_N)
    eps_val = 1e-15

    # Component 1 nodes & log probs
    y1_nodes = mu1 .+ (sqrt(2.0) * sigma1) .* x
    p1_nodes = clamp.(normcdf.(y1_nodes), eps_val, 1.0 - eps_val)
    log_p1 = log.(p1_nodes)
    log_1mp1 = log1p.(-p1_nodes)

    # Component 2 nodes & log probs
    y2_nodes = mu2 .+ (sqrt(2.0) * sigma2) .* x
    p2_nodes = clamp.(normcdf.(y2_nodes), eps_val, 1.0 - eps_val)
    log_p2 = log.(p2_nodes)
    log_1mp2 = log1p.(-p2_nodes)

    log_pi = log(pi_val)
    log_1mpi = log1p(-pi_val)

    K = length(x)
    node_terms1 = Vector{Float64}(undef, K)
    node_terms2 = Vector{Float64}(undef, K)

    tot_nll = 0.0
    N = length(n_arr)

    for i in 1:N
        n_i = n_arr[i]
        h_i = h_arr[i]
        cnt_i = c_arr[i]

        log_c = log_comb(n_i, h_i)

        @inbounds for k in 1:K
            node_terms1[k] = log_w[k] + log_c + h_i * log_p1[k] + (n_i - h_i) * log_1mp1[k]
            node_terms2[k] = log_w[k] + log_c + h_i * log_p2[k] + (n_i - h_i) * log_1mp2[k]
        end

        log_I1 = logsumexp(node_terms1)
        log_I2 = logsumexp(node_terms2)

        log_prob = logaddexp(log_pi + log_I1, log_1mpi + log_I2)
        tot_nll -= cnt_i * log_prob
    end

    return tot_nll
end

"""
    unconstrained_to_mixture_params(params::Vector{Float64})::Tuple{Float64, Float64, Float64, Float64, Float64}

Convert unconstrained vector [theta_mu1, theta_dmu, theta_sig1, theta_sig2, theta_pi]
to physical parameters (mu1, mu2, sigma1, sigma2, pi) with mu1 < mu2.
"""
function unconstrained_to_mixture_params(params::Vector{Float64})::Tuple{Float64, Float64, Float64, Float64, Float64}
    mu1 = params[1]
    dmu = softplus(params[2])
    mu2 = mu1 + dmu
    sigma1 = softplus(params[3])
    sigma2 = softplus(params[4])
    pi_val = logistic(params[5])
    return (mu1, mu2, sigma1, sigma2, pi_val)
end

"""
    mixture_to_unconstrained_params(mu1::Float64, mu2::Float64, sigma1::Float64, sigma2::Float64, pi_val::Float64)::Vector{Float64}

Convert physical parameters (mu1, mu2, sigma1, sigma2, pi) with mu1 < mu2
to unconstrained vector [theta_mu1, theta_dmu, theta_sig1, theta_sig2, theta_pi].
"""
function mixture_to_unconstrained_params(mu1::Float64, mu2::Float64, sigma1::Float64, sigma2::Float64, pi_val::Float64)::Vector{Float64}
    if mu2 <= mu1
        throw(ArgumentError("Require mu1 < mu2 for unconstrained conversion"))
    end
    theta1 = mu1
    theta2 = inv_softplus(mu2 - mu1)
    theta3 = inv_softplus(sigma1)
    theta4 = inv_softplus(sigma2)
    theta5 = logit(pi_val)
    return [theta1, theta2, theta3, theta4, theta5]
end

"""
    nll_probit_normal_mixture_unconstrained(params::Vector{Float64}, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64

Unconstrained parameterization wrapper for Optim.jl.
"""
function nll_probit_normal_mixture_unconstrained(params::Vector{Float64}, n_arr::Vector{Int}, h_arr::Vector{Int}, c_arr::Vector{Int}; GH_N::Int=100)::Float64
    mu1, mu2, sigma1, sigma2, pi_val = unconstrained_to_mixture_params(params)
    return nll_probit_normal_mixture_physical(mu1, mu2, sigma1, sigma2, pi_val, n_arr, h_arr, c_arr; GH_N=GH_N)
end

end # module GaussianMixture
