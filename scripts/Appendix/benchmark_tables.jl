# Shared command-line and text-table helpers for the three rebuttal benchmarks.
# No plotting packages, numerical solver or simulation is loaded here.
using DataFrames, JLD2, Printf, Statistics, LinearAlgebra

function benchmark_options(args, experiment)
    simulate, force, validation_only = "--simulate" in args, "--force" in args, "--validate" in args
    data, tag = "", ""
    for arg in args
        if startswith(arg,"--data=")
            data = split(arg,"=";limit=2)[2]
            isempty(data) && error("--data requires a folder")
        elseif startswith(arg,"--tag=")
            tag = split(arg,"=";limit=2)[2]
            occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$",tag) || error("Invalid run tag")
        elseif !(arg in ["--simulate","--force","--validate","--help"])
            error("Unknown argument: $arg")
        end
    end
    simulate && !isempty(data) && error("--data and --simulate cannot be combined")
    validation_only && !simulate && error("--validate requires --simulate")
    isempty(tag) && (tag = (validation_only ? "validation_" : "run_") * Dates.format(now(),"yyyymmdd_HHMMSSsss"))
    data_dir = isempty(data) ? joinpath(PROJECT_ROOT,"data","published","Appendix",experiment) :
        (isabspath(data) ? data : joinpath(PROJECT_ROOT,data))
    simulate && (data_dir = joinpath(PROJECT_ROOT,"data","generated","Appendix",experiment,tag))
    return (;simulate,force,validation_only,data_dir=normpath(data_dir))
end

function load_benchmark(data_dir, filename)
    file = joinpath(data_dir,filename)
    isfile(file) || error("Saved benchmark missing: $file. Use --simulate explicitly to recompute it.")
    return JLD2.load(file)
end

function print_benchmark_table(io, df, columns)
    println(io,"| ",join(first.(columns)," | ")," |")
    println(io,"| ",join(fill("---",length(columns))," | ")," |")
    for row in eachrow(df)
        values = [row[last(column)] for column in columns]
        println(io,"| ",join([v isa AbstractFloat ? @sprintf("%.12g",v) : string(v) for v in values]," | ")," |")
    end
end

benchmark_run_label(data) = get(data,"validation_only",false) ? "REDUCED VALIDATION" :
    get(data,"generated_by_current_core",false) ? "new run with publication settings" : "archived publication settings"

function print_random_benchmark(data;io=stdout)
    df = data["df_results"]
    println(io,"Under-observed random tasks — ",benchmark_run_label(data))
    print_benchmark_table(io,df,["m"=>:m,"Seed"=>:seed,"TOD"=>:cost_TOD,"LCO"=>:cost_LCO,
        "GD"=>:cost_GD,"LCO improvement (%)"=>:relative_LCO_improvement_over_TOD_percent,
        "LCO–GD gap (%)"=>:relative_LCO_GD_cost_gap_percent])
    println(io,"\nStored run timings (seconds; machine/compiler dependent):")
    timing = copy(df)
    timing.speedup = df.runtime_GD_s ./ df.runtime_LCO_s
    print_benchmark_table(io,timing,["m"=>:m,"TOD (s)"=>:runtime_TOD_s,"LCO (s)"=>:runtime_LCO_s,
        "GD (s)"=>:runtime_GD_s,"GD/LCO time"=>:speedup])
    println(io,"\nDiagnostics (converged denotes the stored final-update residual test):")
    print_benchmark_table(io,df,["m"=>:m,"O rank"=>:observability_rank,
        "LCO passes"=>:LCO_passes,"Final residual"=>:LCO_convergence_residual,
        "LCO converged"=>:LCO_converged,"Max positive increment"=>:LCO_maximum_positive_cost_increment,
        "LCO monotonic"=>:LCO_monotonic,"rho"=>:LCO_selected_rho,
        "GD iterations"=>:GD_iterations,"GD converged"=>:GD_converged,
        "GD hit limit"=>:GD_hit_iteration_limit,"Completed"=>:completed])
    println(io,"\nGD limit flags are retained; cost agreement is not a convergence certificate.")
    println(io,"The first un-warmed LCO call includes possible Julia compilation overhead; the archive does not separate it.")
end

function print_arm_benchmark(data;io=stdout)
    println(io,"Under-observed redundant-arm stabilization — ",benchmark_run_label(data))
    df = data["df_results"]
    print_benchmark_table(io,df,["sigma_eta"=>:sigma_eta,"TOD"=>:cost_TOD,"LCO"=>:cost_LCO,
        "LCO improvement (%)"=>:relative_LCO_improvement_over_TOD_percent,
        "LCO passes"=>:LCO_passes,"Final residual"=>:LCO_convergence_residual,
        "Monotonic"=>:LCO_monotonic,"rho"=>:LCO_selected_rho,"Completed"=>:completed])
    println(io,"\nStored run timings (seconds; machine/compiler dependent):")
    print_benchmark_table(io,df,["sigma_eta"=>:sigma_eta,"TOD (s)"=>:runtime_TOD_s,"LCO (s)"=>:runtime_LCO_s])
    println(io,"\nObservability rank: ",data["observability"].numerical_rank,"/6; singular values: ",data["observability"].singular_values)
    println(io,"No GD or Monte Carlo comparison was performed.")
end

function initialization_summary(data)
    df = data["df_results"]
    reference = only(df.final_analytical_cost[df.initialization_family .== "TOD"])
    cost_TOD = data["study_metadata"].cost_TOD
    families = ["TOD","zero","random_unit_scale","random_large_scale","all_random"]
    rows = NamedTuple[]
    for family in families
        selected = family == "all_random" ? startswith.(df.initialization_family,"random_") : df.initialization_family .== family
        group = df[selected,:]
        final = group.final_analytical_cost
        gaps = 100 .* (final .- reference) ./ reference
        finite = isfinite.(final)
        push!(rows,(;family,n=nrow(group),finite=count(finite),below_TOD=count(final .< cost_TOD),
            min_cost=minimum(final),median_cost=median(final),max_cost=maximum(final),
            min_improvement=minimum(100 .* (cost_TOD .- final) ./ cost_TOD),
            max_improvement=maximum(100 .* (cost_TOD .- final) ./ cost_TOD),
            within_0p1_pct=100*count(abs.(gaps) .<= 0.1)/length(gaps),
            within_0p2_pct=100*count(abs.(gaps) .<= 0.2)/length(gaps),
            max_abs_gap_pct=maximum(abs.(gaps)),
            max_initial_ratio=maximum(group.initial_analytical_cost ./ cost_TOD),
            exceptions=count(occursin.("exception",group.failure_reason)),
            tolerance_reached=count(group.tolerance_reached)))
    end
    return DataFrame(rows)
end

function print_initialization_benchmark(data;io=stdout)
    println(io,"LCO initialization sensitivity — ",benchmark_run_label(data))
    println(io,"Task seed: ",data["study_metadata"].task_seed,"; TOD cost: ",data["study_metadata"].cost_TOD)
    table = initialization_summary(data)
    print_benchmark_table(io,table,["Family"=>:family,"n"=>:n,"Min final cost"=>:min_cost,
        "Median"=>:median_cost,"Max"=>:max_cost,"Below TOD"=>:below_TOD,
        "Within 0.1% (%)"=>:within_0p1_pct,"Within 0.2% (%)"=>:within_0p2_pct,
        "Max abs gap (%)"=>:max_abs_gap_pct,"Max initial/TOD"=>:max_initial_ratio])
    println(io,"\nFinite final costs, solver exceptions and early-stopping diagnostics:")
    print_benchmark_table(io,table,["Family"=>:family,"Finite"=>:finite,"Exceptions"=>:exceptions,
        "Tolerance reached"=>:tolerance_reached])
    println(io,"\nGaps use the TOD-initialized LCO cost; within-threshold counts use absolute gaps.")
    println(io,"Finite decrease and reaching the 1e-10 stopping tolerance are separate diagnostics.")
    println(io,"These results provide empirical robustness evidence, not a global-optimality guarantee.")
end

function observability_diagnostics(A,H,T)
    m = size(A,1)
    observability_matrix = zeros(T*size(H,1),m)
    A_power = Matrix{Float64}(I,m,m)
    for t in 1:T
        observability_matrix[(t-1)*size(H,1)+1:t*size(H,1),:] .= H*A_power
        A_power = A_power*A
    end
    singular_values = svdvals(observability_matrix)
    rank_tolerance = maximum(size(observability_matrix))*eps(Float64)*maximum(singular_values)
    numerical_rank = count(singular_values .> rank_tolerance)
    nonzero = singular_values[singular_values .> rank_tolerance]
    return (;numerical_rank,condition_number=numerical_rank == m ? maximum(singular_values)/minimum(singular_values) : Inf,
        minimum_singular_value=minimum(singular_values),rank_tolerance,singular_values,
        nonzero_condition_number=isempty(nonzero) ? Inf : maximum(nonzero)/minimum(nonzero))
end
