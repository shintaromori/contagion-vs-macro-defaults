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

function run_fixed_parameter_validation()
    data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
    ref_path = joinpath(@__DIR__, "..", "test", "reference_python_results.csv")
    
    if !isfile(data_path)
        error("Data file not found at $data_path")
    end
    if !isfile(ref_path)
        error("Reference file not found at $ref_path")
    end

    df_data = CSV.read(data_path, DataFrame)
    df_ref = CSV.read(ref_path, DataFrame)

    println("==========================================================================")
    println("PHASE J1: FIXED-PARAMETER CROSS-LANGUAGE VALIDATION (JULIA vs PYTHON)")
    println("==========================================================================")
    
    all_pass = true
    max_delta_nll = 0.0

    for row in eachrow(df_ref)
        cls = row.class
        pt_type = row.point_type
        mu = Float64(row.mu)
        sigma = Float64(row.sigma)
        p0 = Float64(row.p0)
        rhoA = Float64(row.rhoA)
        py_nll = Float64(row.python_nll)

        n_col = Symbol(cls)
        h_col = Symbol("D_$cls")

        n_vec = Vector{Int}(df_data[!, n_col])
        h_vec = Vector{Int}(df_data[!, h_col])

        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        # 1. Probit-Normal NLL in Julia
        jl_nll_pn = nll_probit_normal(mu, sigma, n_arr, h_arr, c_arr; GH_N=100)
        delta_pn = abs(jl_nll_pn - py_nll)

        # 2. Vasicek NLL in Julia (via direct Vasicek PMF)
        jl_nll_vas = nll_vasicek(p0, rhoA, n_arr, h_arr, c_arr; GH_N=100)
        delta_vas = abs(jl_nll_vas - py_nll)

        # 3. Equivalence in Julia (mapped Vasicek -> ProbitNormal NLL)
        mu_mapped, sigma_mapped = vasicek_to_probit_normal(p0, rhoA)
        jl_nll_mapped = nll_probit_normal(mu_mapped, sigma_mapped, n_arr, h_arr, c_arr; GH_N=100)
        delta_eq = abs(jl_nll_mapped - jl_nll_pn)

        println("Class: $cls | Point: $pt_type")
        println("  Python Ref NLL : $(rpad(py_nll, 18))")
        println("  Julia Probit NLL: $(rpad(jl_nll_pn, 18)) (Delta vs Python: $(rpad(delta_pn, 12)))")
        println("  Julia Vasicek NLL: $(rpad(jl_nll_vas, 18)) (Delta vs Python: $(rpad(delta_vas, 12)))")
        println("  Julia Equiv NLL  : $(rpad(jl_nll_mapped, 18)) (Vasicek->Probit Delta: $(rpad(delta_eq, 12)))")

        if delta_pn > max_delta_nll
            max_delta_nll = delta_pn
        end

        # Hard numerical gate threshold <= 1e-6
        if delta_pn > 1e-6
            println("  --> [FAIL] Delta NLL exceeds hard gate threshold 1e-6!")
            all_pass = false
        else
            println("  --> [PASS]")
        end
        println("--------------------------------------------------------------------------")
    end

    println("==========================================================================")
    println("SUMMARY: Fixed-Parameter Cross-Language Gate")
    println("  Max Delta NLL across all test cases: $max_delta_nll")
    if all_pass
        println("  STATUS: PASSED (All cases within 1e-6 threshold)")
    else
        println("  STATUS: FAILED (One or more cases exceeded 1e-6 threshold)")
    end
    println("==========================================================================")

    return all_pass
end

if abspath(PROGRAM_FILE) == @__FILE__
    pass = run_fixed_parameter_validation()
    exit(pass ? 0 : 1)
end
