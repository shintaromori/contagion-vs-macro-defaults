using Test
using CSV
using DataFrames
using StatsFuns: normcdf

@testset "5. Phase J2 Gaussian Mixture (M3) Unit & Equivalence Tests" begin

    data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
    @test isfile(data_path)
    df_data = CSV.read(data_path, DataFrame)

    m0_params = Dict(
        "ALL" => (mu = -2.4355093427774697, sigma = 0.3960412471843063, py_m0_nll = 431.3637105779871),
        "SG"  => (mu = -2.033719595595432,  sigma = 0.4194127450554078, py_m0_nll = 417.7085349447972),
        "IG"  => (mu = -3.292821080402438,  sigma = 0.5000005519247075, py_m0_nll = 182.38168604609683)
    )

    m3_best_candidates = Dict(
        "ALL" => (mu1 = -2.580052602925847,  mu2 = -2.4117600823410132, sig1 = 0.6970853207813303, sig2 = 0.376632390070718,   pi = 0.19983508152655866, py_nll = 430.271292090558),
        "SG"  => (mu1 = -3.1958753233822024, mu2 = -1.9688028557365325, sig1 = 1.0793210510643088, sig2 = 0.3424642477500067,  pi = 0.12496466030350036, py_nll = 416.2635553361757),
        "IG"  => (mu1 = -3.4194284031589737, mu2 = -2.608761609066608,  sig1 = 0.4995469023836386, sig2 = 5.1401023355998315e-08, pi = 0.8835765666884274, py_nll = 180.1129479180354)
    )

    @testset "5A. Collapsed Boundary M3 -> M0 & Pi Invariance" begin
        for cls in ["ALL", "SG", "IG"]
            n_vec = Vector{Int}(df_data[!, Symbol(cls)])
            h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
            n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

            mu0 = m0_params[cls].mu
            sig0 = m0_params[cls].sigma

            jl_nll_m0 = nll_probit_normal(mu0, sig0, n_arr, h_arr, c_arr; GH_N=100)
            jl_nll_m3_col5 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.5, n_arr, h_arr, c_arr; GH_N=100)
            jl_nll_m3_col2 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.2, n_arr, h_arr, c_arr; GH_N=100)
            jl_nll_m3_col8 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.8, n_arr, h_arr, c_arr; GH_N=100)

            @test abs(jl_nll_m3_col5 - jl_nll_m0) < 1e-6
            @test abs(jl_nll_m3_col2 - jl_nll_m3_col5) < 1e-6
            @test abs(jl_nll_m3_col8 - jl_nll_m3_col5) < 1e-6
        end
    end

    @testset "5B. Label-Swap Invariance" begin
        for cls in ["ALL", "SG", "IG"]
            n_vec = Vector{Int}(df_data[!, Symbol(cls)])
            h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
            n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

            cand = m3_best_candidates[cls]
            nll_orig = nll_probit_normal_mixture_physical(cand.mu1, cand.mu2, cand.sig1, cand.sig2, cand.pi, n_arr, h_arr, c_arr; GH_N=100)
            nll_swap = nll_probit_normal_mixture_physical(cand.mu2, cand.mu1, cand.sig2, cand.sig1, 1.0 - cand.pi, n_arr, h_arr, c_arr; GH_N=100)

            @test abs(nll_orig - nll_swap) < 1e-6
        end
    end

    @testset "5C. Python vs Julia Fixed-Parameter M3 Likelihood" begin
        for cls in ["ALL", "SG", "IG"]
            n_vec = Vector{Int}(df_data[!, Symbol(cls)])
            h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
            n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

            cand = m3_best_candidates[cls]
            jl_nll = nll_probit_normal_mixture_physical(cand.mu1, cand.mu2, cand.sig1, cand.sig2, cand.pi, n_arr, h_arr, c_arr; GH_N=100)
            @test abs(jl_nll - cand.py_nll) < 1e-6
        end
    end

    @testset "5D. IG Near-Degenerate Sigma (5.14e-8) Stability" begin
        n_vec_ig = Vector{Int}(df_data[!, :IG])
        h_vec_ig = Vector{Int}(df_data[!, :D_IG])
        n_arr_ig, h_arr_ig, c_arr_ig = group_nh(n_vec_ig, h_vec_ig)

        ig_cand = m3_best_candidates["IG"]
        jl_ig_nll = nll_probit_normal_mixture_physical(
            ig_cand.mu1, ig_cand.mu2, ig_cand.sig1, ig_cand.sig2, ig_cand.pi,
            n_arr_ig, h_arr_ig, c_arr_ig; GH_N=100
        )
        @test isfinite(jl_ig_nll)
        @test abs(jl_ig_nll - ig_cand.py_nll) < 1e-6
    end

    @testset "5E. Unconstrained Parameter Transformation Round-Trip" begin
        params_unconstrained = [-2.5, 0.5, -0.4, -1.0, 0.0]
        mu1, mu2, sig1, sig2, pi_val = unconstrained_to_mixture_params(params_unconstrained)
        @test mu1 < mu2
        @test sig1 > 0.0
        @test sig2 > 0.0
        @test 0.0 < pi_val < 1.0
    end

    @testset "5F. Physical <-> Unconstrained Round-Trip Tests (ALL, SG, IG)" begin
        for cls in ["ALL", "SG", "IG"]
            cand = m3_best_candidates[cls]
            x0 = mixture_to_unconstrained_params(cand.mu1, cand.mu2, cand.sig1, cand.sig2, cand.pi)
            mu1_rec, mu2_rec, sig1_rec, sig2_rec, pi_rec = unconstrained_to_mixture_params(x0)

            @test isapprox(cand.mu1, mu1_rec, atol=1e-10)
            @test isapprox(cand.mu2, mu2_rec, atol=1e-10)
            @test isapprox(cand.sig1, sig1_rec, atol=1e-10)
            @test isapprox(cand.sig2, sig2_rec, atol=1e-10)
            @test isapprox(cand.pi, pi_rec, atol=1e-10)
        end
    end

end
