using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "IDMLikelihoods.jl"))
using .IDMLikelihoods

using CSV
using DataFrames
using Test

function group_nh(n_series::Vector{Int}, h_series::Vector{Int})
    pairs = Dict{Tuple{Int, Int}, Int}()
    for (n, h) in zip(n_series, h_series)
        pairs[(n, h)] = get(pairs, (n, h), 0) + 1
    end
    n_arr = Int[]
    h_arr = Int[]
    c_arr = Int[]
    for ((n, h), cnt) in pairs
        push!(n_arr, n)
        push!(h_arr, h)
        push!(c_arr, cnt)
    end
    return n_arr, h_arr, c_arr
end

function run_mixture_fixed_parameter_validation()
    data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
    if !isfile(data_path)
        error("Data file not found at $data_path")
    end

    df_data = CSV.read(data_path, DataFrame)

    println("==========================================================================")
    println("PHASE J2: GAUSSIAN MIXTURE (M3) FIXED-PARAMETER VALIDATION")
    println("==========================================================================")

    all_pass = true

    # M0 reference parameters from validated J1/Phase 2A results
    m0_params = Dict(
        "ALL" => (mu = -2.4355093427774697, sigma = 0.3960412471843063, py_m0_nll = 431.3637105779871),
        "SG"  => (mu = -2.033719595595432,  sigma = 0.4194127450554078, py_m0_nll = 417.7085349447972),
        "IG"  => (mu = -3.292821080402438,  sigma = 0.5000005519247075, py_m0_nll = 182.38168604609683)
    )

    # 1. Collapsed Boundary & Pi Invariance Test
    println("\n[1] Collapsed Boundary & Pi Invariance Tests (Physical Likelihood)")
    println("--------------------------------------------------------------------------")

    collapsed_results = Dict()

    for cls in ["ALL", "SG", "IG"]
        n_vec = Vector{Int}(df_data[!, Symbol(cls)])
        h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        mu0 = m0_params[cls].mu
        sig0 = m0_params[cls].sigma
        py_nll0 = m0_params[cls].py_m0_nll

        # Evaluate M0
        jl_nll_m0 = nll_probit_normal(mu0, sig0, n_arr, h_arr, c_arr; GH_N=100)

        # Evaluate Collapsed M3 at pi = 0.5
        jl_nll_m3_col5 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.5, n_arr, h_arr, c_arr; GH_N=100)

        # Evaluate Collapsed M3 at pi = 0.2 and pi = 0.8
        jl_nll_m3_col2 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.2, n_arr, h_arr, c_arr; GH_N=100)
        jl_nll_m3_col8 = nll_probit_normal_mixture_physical(mu0, mu0, sig0, sig0, 0.8, n_arr, h_arr, c_arr; GH_N=100)

        delta_col = abs(jl_nll_m3_col5 - jl_nll_m0)
        delta_pi2 = abs(jl_nll_m3_col2 - jl_nll_m3_col5)
        delta_pi8 = abs(jl_nll_m3_col8 - jl_nll_m3_col5)

        collapsed_results[cls] = (delta_col = delta_col, delta_pi2 = delta_pi2, delta_pi8 = delta_pi8)

        println("Class: $cls")
        println("  M0 NLL           : $jl_nll_m0")
        println("  Collapsed M3 (pi=0.5): $jl_nll_m3_col5 (Delta vs M0: $delta_col)")
        println("  Collapsed M3 (pi=0.2): $jl_nll_m3_col2 (Delta vs pi=0.5: $delta_pi2)")
        println("  Collapsed M3 (pi=0.8): $jl_nll_m3_col8 (Delta vs pi=0.5: $delta_pi8)")

        pass_col = (delta_col <= 1e-6) && (delta_pi2 <= 1e-6) && (delta_pi8 <= 1e-6)
        println("  --> STATUS: $(pass_col ? "[PASS]" : "[FAIL]")")
        if !pass_col
            all_pass = false
        end
    end

    # 2. Label-Swap Invariance Test
    println("\n[2] Label-Swap Invariance Tests (Physical Likelihood)")
    println("--------------------------------------------------------------------------")

    # M3 Best candidates from Python Phase 2B CSVs
    m3_best_candidates = Dict(
        "ALL" => (mu1 = -2.580052602925847,  mu2 = -2.4117600823410132, sig1 = 0.6970853207813303, sig2 = 0.376632390070718,   pi = 0.19983508152655866, py_nll = 430.271292090558),
        "SG"  => (mu1 = -3.1958753233822024, mu2 = -1.9688028557365325, sig1 = 1.0793210510643088, sig2 = 0.3424642477500067,  pi = 0.12496466030350036, py_nll = 416.2635553361757),
        "IG"  => (mu1 = -3.4194284031589737, mu2 = -2.608761609066608,  sig1 = 0.4995469023836386, sig2 = 5.1401023355998315e-08, pi = 0.8835765666884274, py_nll = 180.1129479180354)
    )

    label_swap_results = Dict()

    for cls in ["ALL", "SG", "IG"]
        n_vec = Vector{Int}(df_data[!, Symbol(cls)])
        h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        cand = m3_best_candidates[cls]
        mu1, mu2, sig1, sig2, pi_val = cand.mu1, cand.mu2, cand.sig1, cand.sig2, cand.pi

        # Original order
        nll_orig = nll_probit_normal_mixture_physical(mu1, mu2, sig1, sig2, pi_val, n_arr, h_arr, c_arr; GH_N=100)

        # Swapped order: (mu2, sig2, mu1, sig1, 1 - pi)
        nll_swap = nll_probit_normal_mixture_physical(mu2, mu1, sig2, sig1, 1.0 - pi_val, n_arr, h_arr, c_arr; GH_N=100)

        delta_swap = abs(nll_orig - nll_swap)
        label_swap_results[cls] = delta_swap

        println("Class: $cls")
        println("  Original Vector NLL : $nll_orig")
        println("  Swapped Vector NLL  : $nll_swap (Delta: $delta_swap)")
        pass_swap = delta_swap <= 1e-6
        println("  --> STATUS: $(pass_swap ? "[PASS]" : "[FAIL]")")
        if !pass_swap
            all_pass = false
        end
    end

    # 3. Fixed-Parameter Python vs Julia M3 Validation
    println("\n[3] Python vs Julia Fixed-Parameter M3 Likelihood Audit")
    println("--------------------------------------------------------------------------")

    m3_delta_results = Dict()

    for cls in ["ALL", "SG", "IG"]
        n_vec = Vector{Int}(df_data[!, Symbol(cls)])
        h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        cand = m3_best_candidates[cls]
        jl_nll = nll_probit_normal_mixture_physical(cand.mu1, cand.mu2, cand.sig1, cand.sig2, cand.pi, n_arr, h_arr, c_arr; GH_N=100)
        delta_py = abs(jl_nll - cand.py_nll)
        m3_delta_results[cls] = delta_py

        println("Class: $cls | Best-Found M3 Candidate")
        println("  Python Ref NLL : $(cand.py_nll)")
        println("  Julia M3 NLL   : $jl_nll (Delta: $delta_py)")
        pass_py = delta_py <= 1e-6
        println("  --> STATUS: $(pass_py ? "[PASS]" : "[FAIL]")")
        if !pass_py
            all_pass = false
        end
    end

    # 4. IG Near-Degenerate Candidate Stability Test
    println("\n[4] IG Near-Degenerate Candidate (sigma2 ≈ 5.14e-8) Stability Check")
    println("--------------------------------------------------------------------------")
    ig_cand = m3_best_candidates["IG"]
    n_vec_ig = Vector{Int}(df_data[!, :IG])
    h_vec_ig = Vector{Int}(df_data[!, :D_IG])
    n_arr_ig, h_arr_ig, c_arr_ig = group_nh(n_vec_ig, h_vec_ig)

    jl_ig_nll = nll_probit_normal_mixture_physical(
        ig_cand.mu1, ig_cand.mu2, ig_cand.sig1, ig_cand.sig2, ig_cand.pi,
        n_arr_ig, h_arr_ig, c_arr_ig; GH_N=100
    )
    ig_is_finite = isfinite(jl_ig_nll)
    ig_delta = abs(jl_ig_nll - ig_cand.py_nll)

    println("IG Near-Degenerate NLL: $jl_ig_nll (is_finite: $ig_is_finite, Delta vs Python: $ig_delta)")
    pass_ig = ig_is_finite && (ig_delta <= 1e-6)
    println("  --> STATUS: $(pass_ig ? "[PASS]" : "[FAIL]")")
    if !pass_ig
        all_pass = false
    end

    println("==========================================================================")
    println("SUMMARY: Phase J2 Fixed-Parameter & Equivalence Gate")
    if all_pass
        println("  FINAL STATUS: PASSED (All tests within 1e-6 threshold)")
    else
        println("  FINAL STATUS: FAILED")
    end
    println("==========================================================================")

    return all_pass
end

if abspath(PROGRAM_FILE) == @__FILE__
    pass = run_mixture_fixed_parameter_validation()
    exit(pass ? 0 : 1)
end
