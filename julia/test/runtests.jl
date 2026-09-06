using Test
using CSV
using DataFrames
using StatsFuns: normcdf, norminvcdf

include(joinpath(@__DIR__, "..", "src", "IDMLikelihoods.jl"))
using .IDMLikelihoods

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

@testset "IDMLikelihoods Test Suite (Phase J1 & Phase J2)" begin

    @testset "1. Hermite Quadrature Utilities" begin
        x, log_w = get_hermite_quadrature(100)
        @test length(x) == 100
        @test length(log_w) == 100
        # Weights must sum to 1.0 (normalized)
        @test isapprox(sum(exp.(log_w)), 1.0, atol=1e-12)
        # Symmetry of nodes
        @test isapprox(x[1], -x[end], atol=1e-12)
    end

    @testset "2. Log Binomial Coefficient" begin
        @test isapprox(log_comb(10, 0), 0.0, atol=1e-12)
        @test isapprox(log_comb(10, 10), 0.0, atol=1e-12)
        @test isapprox(log_comb(5, 2), log(10.0), atol=1e-12)
    end

    @testset "3. Parameter Round-Trip Transformations" begin
        mu_in, sigma_in = -2.435509, 0.3960413
        p0, rhoA = probit_normal_to_vasicek(mu_in, sigma_in)
        mu_out, sigma_out = vasicek_to_probit_normal(p0, rhoA)
        @test isapprox(mu_in, mu_out, atol=1e-12)
        @test isapprox(sigma_in, sigma_out, atol=1e-12)
    end

    @testset "4. Cross-Language Fixed-Parameter Validation" begin
        data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
        ref_path = joinpath(@__DIR__, "reference_python_results.csv")

        @test isfile(data_path)
        @test isfile(ref_path)

        df_data = CSV.read(data_path, DataFrame)
        df_ref = CSV.read(ref_path, DataFrame)

        for row in eachrow(df_ref)
            cls = row.class
            mu = Float64(row.mu)
            sigma = Float64(row.sigma)
            py_nll = Float64(row.python_nll)

            n_col = Symbol(cls)
            h_col = Symbol("D_$cls")

            n_vec = Vector{Int}(df_data[!, n_col])
            h_vec = Vector{Int}(df_data[!, h_col])

            n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

            # Check Probit-Normal NLL match vs Python reference
            jl_nll_pn = nll_probit_normal(mu, sigma, n_arr, h_arr, c_arr; GH_N=100)
            @test abs(jl_nll_pn - py_nll) < 1e-6
        end
    end

    @testset "5. Julia-Internal Direct Vasicek vs Probit-Normal Likelihood Equivalence" begin
        data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
        ref_path = joinpath(@__DIR__, "reference_python_results.csv")

        df_data = CSV.read(data_path, DataFrame)
        df_ref = CSV.read(ref_path, DataFrame)

        for row in eachrow(df_ref[df_ref.point_type .== "MLE", :])
            cls = row.class
            mu = Float64(row.mu)
            sigma = Float64(row.sigma)
            p0 = Float64(row.p0)
            rhoA = Float64(row.rhoA)

            n_col = Symbol(cls)
            h_col = Symbol("D_$cls")

            n_vec = Vector{Int}(df_data[!, n_col])
            h_vec = Vector{Int}(df_data[!, h_col])

            n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

            jl_nll_pn = nll_probit_normal(mu, sigma, n_arr, h_arr, c_arr; GH_N=100)
            jl_nll_vas = nll_vasicek(p0, rhoA, n_arr, h_arr, c_arr; GH_N=100)

            delta_direct = abs(jl_nll_vas - jl_nll_pn)
            @test delta_direct < 1e-6
        end
    end

    include("test_gaussian_mixture.jl")

end
