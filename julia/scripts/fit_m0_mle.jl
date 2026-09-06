using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "IDMLikelihoods.jl"))
using .IDMLikelihoods

using CSV
using DataFrames
using Dates
using Optim
using StatsFuns: norminvcdf

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

function fit_m0_all_classes()
    data_path = joinpath(@__DIR__, "..", "..", "data", "M_1920_2023.csv")
    ref_path = joinpath(@__DIR__, "..", "test", "reference_python_results.csv")

    df_data = CSV.read(data_path, DataFrame)
    df_ref = CSV.read(ref_path, DataFrame)

    # Reference Python MLE values
    py_mle_targets = Dict(
        "ALL" => (nll=431.3637105785, mu=-2.435509, sigma=0.3960413),
        "SG"  => (nll=417.7085349450, mu=-2.033720, sigma=0.4194128),
        "IG"  => (nll=182.3816860461, mu=-3.292821, sigma=0.5000005)
    )

    results_rows = []
    benchmark_rows = []

    println("==========================================================================")
    println("PHASE J1: M0 MLE OPTIMIZATION & COMPARISON (JULIA vs PYTHON)")
    println("==========================================================================")

    all_pass = true

    for cls in ["ALL", "SG", "IG"]
        n_col = Symbol(cls)
        h_col = Symbol("D_$cls")

        n_vec = Vector{Int}(df_data[!, n_col])
        h_vec = Vector{Int}(df_data[!, h_col])

        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        # Initial point candidates
        p_hat = sum(h_vec) / sum(n_vec)
        mu_hat = norminvcdf(clamp(p_hat, 1e-5, 1.0 - 1e-5))

        x0_candidates = [
            [mu_hat, log(0.3)],
            [mu_hat, log(0.1)],
            [mu_hat - 0.1, log(0.3)],
            [mu_hat, log(0.4)]
        ]

        # Objective function for Optim.jl
        f_obj(th) = nll_probit_normal_unconstrained(th, n_arr, h_arr, c_arr; GH_N=100)

        # Start benchmark timer
        t_start = now()
        t_start_sec = time()

        best_res = nothing
        best_nll = Inf

        start_count = length(x0_candidates)
        total_iters = 0
        total_evals = 0

        for x0 in x0_candidates
            res = Optim.optimize(f_obj, x0, LBFGS(), Optim.Options(g_tol=1e-8, iterations=10000))
            total_iters += Optim.iterations(res)
            total_evals += Optim.f_calls(res)
            if res.minimum < best_nll
                best_nll = res.minimum
                best_res = res
            end
        end

        t_end_sec = time()
        t_end = now()
        elapsed = t_end_sec - t_start_sec

        mu_fit = best_res.minimizer[1]
        sigma_fit = exp(best_res.minimizer[2])
        jl_nll = best_res.minimum

        p0_eq, rhoA_eq = probit_normal_to_vasicek(mu_fit, sigma_fit)

        py_ref = py_mle_targets[cls]
        delta_nll = abs(jl_nll - py_ref.nll)
        delta_mu = abs(mu_fit - py_ref.mu)
        delta_sigma = abs(sigma_fit - py_ref.sigma)

        is_converged = Optim.converged(best_res)

        println("Class: $cls")
        println("  Python Ref MLE NLL : $(py_ref.nll)")
        println("  Julia MLE NLL      : $(jl_nll) (Delta: $delta_nll)")
        println("  Python Ref Params  : mu = $(py_ref.mu), sigma = $(py_ref.sigma)")
        println("  Julia Params       : mu = $mu_fit, sigma = $sigma_fit (p0_eq = $p0_eq, rhoA_eq = $rhoA_eq)")
        println("  Convergence Status : $(is_converged ? "CONVERGED" : "FAILED") (total_starts: $start_count, total_iters: $total_iters, total_f_calls: $total_evals, time: $(round(elapsed, digits=4))s)")

        pass_cls = (delta_nll <= 1e-6) && is_converged
        if pass_cls
            println("  --> [PASS]")
        else
            println("  --> [FAIL]")
            all_pass = false
        end
        println("--------------------------------------------------------------------------")

        # Record validation result row
        push!(results_rows, (
            class = cls,
            model = "M0 (Probit-Normal NULL)",
            python_ref_nll = py_ref.nll,
            julia_mle_nll = jl_nll,
            delta_nll = delta_nll,
            julia_mu = mu_fit,
            julia_sigma = sigma_fit,
            julia_p0_equiv = p0_eq,
            julia_rhoA_equiv = rhoA_eq,
            converged = is_converged,
            status_pass = pass_cls
        ))

        # Record benchmark row
        push!(benchmark_rows, (
            timestamp = Dates.format(t_start, "yyyy-mm-ddTHH:MM:SS"),
            julia_version = string(VERSION),
            class = cls,
            model = "M0",
            GH_N = 100,
            optimizer = "LBFGS",
            start_count = start_count,
            iterations = total_iters,
            evaluations = total_evals,
            status = is_converged ? "converged" : "failed",
            elapsed_seconds = round(elapsed, digits=6),
            nll = jl_nll,
            mu_hat = mu_fit,
            sigma_hat = sigma_fit
        ))
    end

    println("==========================================================================")
    println("SUMMARY: Julia M0 MLE Optimization & Benchmark Export")
    if all_pass
        println("  STATUS: PASSED (All classes matched Python MLE within 1e-6)")
    else
        println("  STATUS: FAILED")
    end
    println("==========================================================================")

    # Save validation & benchmark outputs into pdata/
    pdata_dir = joinpath(@__DIR__, "..", "..", "pdata")
    val_out_path = joinpath(pdata_dir, "revision2026_julia_j1_validation_results.csv")
    bench_out_path = joinpath(pdata_dir, "revision2026_julia_j1_benchmark.csv")

    CSV.write(val_out_path, DataFrame(results_rows))
    CSV.write(bench_out_path, DataFrame(benchmark_rows))

    println("Saved validation output to: $val_out_path")
    println("Saved benchmark output to: $bench_out_path")

    return all_pass
end

if abspath(PROGRAM_FILE) == @__FILE__
    pass = fit_m0_all_classes()
    exit(pass ? 0 : 1)
end
