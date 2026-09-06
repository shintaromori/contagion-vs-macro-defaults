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

function compute_aic_bic(nll::Float64, k::Int, T::Int)
    logL = -nll
    aic = 2.0 * k - 2.0 * logL
    bic = Float64(k) * log(Float64(T)) - 2.0 * logL
    return (logL = logL, aic = aic, bic = bic)
end

function run_full_phase_j2_reproduction()
    t_wall_start = time()

    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    pdata_dir = joinpath(project_root, "pdata")
    if !isdir(pdata_dir)
        mkpath(pdata_dir)
    end

    data_path = joinpath(project_root, "data", "M_1920_2023.csv")
    if !isfile(data_path)
        error("Data file not found at $data_path")
    end

    df_data = CSV.read(data_path, DataFrame)

    # Reference Python Phase 2B best-found values
    py_m3_targets = Dict(
        "ALL" => (nll=430.271292090558, mu1=-2.580052602925847,  mu2=-2.4117600823410132, sig1=0.6970853207813303, sig2=0.376632390070718,   pi=0.19983508152655866),
        "SG"  => (nll=416.177207454587, mu1=-2.837252756919391, mu2=-1.9626870108228043, sig1=0.8958485309711917, sig2=0.33744589381289813, pi=0.16108743264688935),
        "IG"  => (nll=180.1129479180354, mu1=-3.4194284031589737, mu2=-2.608761609066608,  sig1=0.4995469023836386, sig2=5.1401023355998315e-08, pi=0.8835765666884274)
    )

    m0_params = Dict(
        "ALL" => (mu = -2.4355093427774697, sigma = 0.3960412471843063, nll = 431.3637105780085),
        "SG"  => (mu = -2.033719595595432,  sigma = 0.4194127450554078, nll = 417.708534944726),
        "IG"  => (mu = -3.292821080402438,  sigma = 0.5000005519247075, nll = 182.38168604609112)
    )

    println("==========================================================================")
    println("PHASE J2: FULL 540-START MULTI-START OPTIMIZATION REPRODUCTION")
    println("==========================================================================")

    # ----------------------------------------------------------------------
    # JIT WARM-UP
    # ----------------------------------------------------------------------
    println("\n[0] Explicit JIT Warm-up Run...")
    n_vec_all = Vector{Int}(df_data[!, :ALL])
    h_vec_all = Vector{Int}(df_data[!, :D_ALL])
    n_arr_w, h_arr_w, c_arr_w = group_nh(n_vec_all, h_vec_all)

    t_warmup_start = time()
    f_warmup(th) = nll_probit_normal_mixture_unconstrained(th, n_arr_w, h_arr_w, c_arr_w; GH_N=100)
    x0_warmup = mixture_to_unconstrained_params(-2.5, -2.4, 0.6, 0.4, 0.3)
    res_warmup = Optim.optimize(f_warmup, x0_warmup, NelderMead(), Optim.Options(iterations=5))
    t_warmup_end = time()
    warmup_time = t_warmup_end - t_warmup_start
    println("  JIT Warm-up completed in $(round(warmup_time, digits=4))s (Excluded from benchmarks).")

    diagnostics_rows = []
    best_estimates_rows = []
    model_comparison_rows = []
    benchmark_rows = []

    total_steady_opt_time = 0.0
    total_iters_all_runs = 0
    total_f_calls_all_runs = 0

    class_timing = Dict()

    # ----------------------------------------------------------------------
    # FULL 180 STARTS PER CLASS (540 TOTAL STARTS)
    # ----------------------------------------------------------------------
    for cls in ["ALL", "SG", "IG"]
        println("\n--------------------------------------------------------------------------")
        println("CLASS: $cls | Executing 180 Multi-Starts...")
        println("--------------------------------------------------------------------------")

        n_vec = Vector{Int}(df_data[!, Symbol(cls)])
        h_vec = Vector{Int}(df_data[!, Symbol("D_$cls")])
        n_arr, h_arr, c_arr = group_nh(n_vec, h_vec)

        mu0 = m0_params[cls].mu
        sig0 = m0_params[cls].sigma

        # Build canonical 180-start grid (5 pi * 4 mu_offset * 3 sig1 * 3 sig2)
        pi_list = [0.1, 0.3, 0.5, 0.7, 0.9]
        mu_offsets = [0.1, 0.5, 1.0, 1.5]
        sig_factor1 = [0.7, 1.0, 1.3]
        sig_factor2 = [0.7, 1.0, 1.3]

        start_grid = []
        for pi_i in pi_list
            for offset in mu_offsets
                mu1_i = mu0 - offset / 2.0
                mu2_i = mu0 + offset / 2.0
                for sf1 in sig_factor1
                    for sf2 in sig_factor2
                        sig1_i = max(1e-3, sig0 * sf1)
                        sig2_i = max(1e-3, sig0 * sf2)
                        push!(start_grid, (pi=pi_i, mu1=mu1_i, sig1=sig1_i, mu2=mu2_i, sig2=sig2_i))
                    end
                end
            end
        end

        f_obj(th) = nll_probit_normal_mixture_unconstrained(th, n_arr, h_arr, c_arr; GH_N=100)

        best_nll = Inf
        best_params = nothing
        second_best_nll = Inf
        second_best_params = nothing

        cls_opt_time = 0.0
        cls_iters = 0
        cls_f_calls = 0
        successful_starts = 0

        for (start_id, sp) in enumerate(start_grid)
            x0 = mixture_to_unconstrained_params(sp.mu1, sp.mu2, sp.sig1, sp.sig2, sp.pi)
            init_nll = f_obj(x0)

            t0 = time()
            res = Optim.optimize(f_obj, x0, NelderMead(), Optim.Options(iterations=6000))
            t1 = time()
            elapsed = t1 - t0

            cls_opt_time += elapsed
            final_nll = res.minimum
            final_mu1, final_mu2, final_sig1, final_sig2, final_pi = unconstrained_to_mixture_params(res.minimizer)
            is_conv = Optim.converged(res)
            iters = Optim.iterations(res)
            f_calls = Optim.f_calls(res)

            cls_iters += iters
            cls_f_calls += f_calls

            if is_conv && isfinite(final_nll)
                successful_starts += 1
            end

            # Mode tracking
            if final_nll < best_nll
                if best_nll < Inf && abs(best_nll - final_nll) > 1e-4
                    second_best_nll = best_nll
                    second_best_params = best_params
                end
                best_nll = final_nll
                best_params = (mu1=final_mu1, mu2=final_mu2, sig1=final_sig1, sig2=final_sig2, pi=final_pi)
            elseif final_nll > best_nll + 1e-4 && final_nll < second_best_nll
                second_best_nll = final_nll
                second_best_params = (mu1=final_mu1, mu2=final_mu2, sig1=final_sig1, sig2=final_sig2, pi=final_pi)
            end

            push!(diagnostics_rows, (
                class = cls,
                start_id = start_id,
                initial_pi = sp.pi,
                initial_mu1 = sp.mu1,
                initial_sigma1 = sp.sig1,
                initial_mu2 = sp.mu2,
                initial_sigma2 = sp.sig2,
                initial_nll = init_nll,
                final_pi = final_pi,
                final_mu1 = final_mu1,
                final_sigma1 = final_sig1,
                final_mu2 = final_mu2,
                final_sigma2 = final_sig2,
                final_nll = final_nll,
                converged = is_conv,
                iterations = iters,
                f_calls = f_calls,
                opt_time_sec = round(elapsed, digits=5)
            ))
        end

        total_steady_opt_time += cls_opt_time
        total_iters_all_runs += cls_iters
        total_f_calls_all_runs += cls_f_calls
        class_timing[cls] = cls_opt_time

        py_target = py_m3_targets[cls].nll
        delta_best = best_nll - py_target

        println("Class $cls Summary:")
        println("  Successful Starts : $successful_starts / 180")
        println("  Best Discovered NLL: $best_nll (Py Ref: $py_target, Delta: $delta_best)")
        println("  Best Parameters    : mu1=$(best_params.mu1), mu2=$(best_params.mu2), sig1=$(best_params.sig1), sig2=$(best_params.sig2), pi=$(best_params.pi)")
        if second_best_nll < Inf
            println("  2nd Mode NLL       : $second_best_nll (Delta vs Best: $(second_best_nll - best_nll))")
            println("  2nd Mode Params    : mu1=$(second_best_params.mu1), mu2=$(second_best_params.mu2), sig1=$(second_best_params.sig1), sig2=$(second_best_params.sig2), pi=$(second_best_params.pi)")
        end
        println("  Class Execution Time: $(round(cls_opt_time, digits=2))s (Iters: $cls_iters, f_calls: $cls_f_calls)")

        # Record parameter estimates
        push!(best_estimates_rows, (class=cls, model="M3 (2-comp Mixture NULL)", parameter="mu1", estimate=best_params.mu1))
        push!(best_estimates_rows, (class=cls, model="M3 (2-comp Mixture NULL)", parameter="mu2", estimate=best_params.mu2))
        push!(best_estimates_rows, (class=cls, model="M3 (2-comp Mixture NULL)", parameter="sigma1", estimate=best_params.sig1))
        push!(best_estimates_rows, (class=cls, model="M3 (2-comp Mixture NULL)", parameter="sigma2", estimate=best_params.sig2))
        push!(best_estimates_rows, (class=cls, model="M3 (2-comp Mixture NULL)", parameter="pi", estimate=best_params.pi))

        # Model comparison (M0 vs M3)
        T_val = 104
        m0_nll_val = m0_params[cls].nll
        m0_crit = compute_aic_bic(m0_nll_val, 2, T_val)
        m3_crit = compute_aic_bic(best_nll, 5, T_val)

        push!(model_comparison_rows, (
            class=cls, model="M0 (Probit-Normal NULL)", T=T_val, k=2,
            logL=m0_crit.logL, nll=m0_nll_val, AIC=m0_crit.aic, BIC=m0_crit.bic
        ))
        push!(model_comparison_rows, (
            class=cls, model="M3 (2-comp Mixture NULL)", T=T_val, k=5,
            logL=m3_crit.logL, nll=best_nll, AIC=m3_crit.aic, BIC=m3_crit.bic
        ))

        # Record benchmark
        push!(benchmark_rows, (
            timestamp = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
            julia_version = string(VERSION),
            class = cls,
            model = "M3",
            GH_N = 100,
            optimizer = "NelderMead",
            start_count = 180,
            successful_starts = successful_starts,
            iterations = cls_iters,
            evaluations = cls_f_calls,
            opt_time_sec = round(cls_opt_time, digits=4),
            best_nll = best_nll,
            mu1_hat = best_params.mu1,
            mu2_hat = best_params.mu2,
            sigma1_hat = best_params.sig1,
            sigma2_hat = best_params.sig2,
            pi_hat = best_params.pi
        ))
    end

    t_wall_end = time()
    total_wall_time = t_wall_end - t_wall_start

    println("\n==========================================================================")
    println("PHASE J2 FULL REPRODUCTION RUNTIME SUMMARY")
    println("==========================================================================")
    println("  JIT Warm-up Time              : $(round(warmup_time, digits=4))s")
    println("  ALL 180-Start Optimization Time: $(round(class_timing["ALL"], digits=2))s")
    println("  SG  180-Start Optimization Time: $(round(class_timing["SG"], digits=2))s")
    println("  IG  180-Start Optimization Time: $(round(class_timing["IG"], digits=2))s")
    println("  Total 540-Start Steady Opt Time: $(round(total_steady_opt_time, digits=2))s")
    println("  Total End-to-End Wall-Clock Time: $(round(total_wall_time, digits=2))s")
    println("  Total Optimizer Iterations    : $total_iters_all_runs")
    println("  Total Function Evaluations    : $total_f_calls_all_runs")
    println("==========================================================================")

    # Save dedicated J2 CSV outputs to pdata/
    diag_csv_path = joinpath(pdata_dir, "revision2026_julia_j2_mixture_diagnostics.csv")
    param_csv_path = joinpath(pdata_dir, "revision2026_julia_j2_mixture_parameter_estimates.csv")
    comp_csv_path = joinpath(pdata_dir, "revision2026_julia_j2_mixture_model_comparison.csv")
    bench_csv_path = joinpath(pdata_dir, "revision2026_julia_j2_benchmark.csv")

    CSV.write(diag_csv_path, DataFrame(diagnostics_rows))
    CSV.write(param_csv_path, DataFrame(best_estimates_rows))
    CSV.write(comp_csv_path, DataFrame(model_comparison_rows))
    CSV.write(bench_csv_path, DataFrame(benchmark_rows))

    println("Saved Phase J2 outputs:")
    println("  - Diagnostics : $diag_csv_path")
    println("  - Parameters  : $param_csv_path")
    println("  - Model Comp  : $comp_csv_path")
    println("  - Benchmark   : $bench_csv_path")

    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    pass = run_full_phase_j2_reproduction()
    exit(pass ? 0 : 1)
end
