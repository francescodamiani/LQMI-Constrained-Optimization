# Shared adapters to the validated library; no duplicated stochastic dynamics.
module IdentifiabilitySimulation
using LinearAlgebra, Random, Statistics, JLD2
include(joinpath(@__DIR__,"..","..","..","src","MyFunctionsOptimalRegularizedSolutions.jl"))
include(joinpath(@__DIR__,"..","..","Main","Figure2","fig2_data.jl"))
include(joinpath(@__DIR__,"..","..","Main","Figure2","reaching_task.jl"))
include(joinpath(@__DIR__,"..","..","Main","Figure2","fig2_statistics.jl"))

function sample_identifiability(seed, N, n, nz, p, T, mx, mz, Σx, Σz,
    A, B, C, H, D, K, M, L, Ωξ, Ωω, Ωη, Q, R)
    # Dataset arrays use (dimension, trial, time); the core returns
    # (dimension, time, trial). Each condition resets its own MersenneTwister.
    mc = monte_carlo_simulation(L,K,A,B,H,C,D,T,n,size(B,2),p,mx,mz,
        Σx,Σz,Q,R,Matrix(Ωξ),Matrix(Ωω),Matrix(Ωη),N;rng=MersenneTwister(seed))
    return permutedims(mc.x_traces,(1,3,2)), permutedims(mc.z_traces,(1,3,2)),
        nothing, mc.empirical_cost, mc.empirical_cost_SE
end

function save_results(validation_only, path, args...)
    results = Dict{String,Any}(args[i]=>args[i+1] for i in 1:2:length(args))
    merge!(results,Dict("validation_only"=>validation_only,"generated_by_current_core"=>true,
        "N_iter_Todorov"=>(validation_only ? 3 : 500),"N_iter_tot"=>(validation_only ? 3 : 500),
        "rho_list"=>[0.0],"rho_list_large"=>[1e-10,1e-8,1e-6],
        "iters_K"=>1,"iters_L"=>1,"use_GD_refinement"=>false,"use_GD_from_IC"=>false))
    JLD2.jldopen(path,"w";compress=true) do file
        for (key,value) in results; file[key]=value; end
    end
    return results
end

function validation_reference()
    # Calibrated reaching dynamics, costs and external-noise amplitudes.
    # Only horizon, trials, iterations and sweep sizes are reduced for validation.
    p = reaching_problem(T_move=0.07)
    d = Dict{String,Any}(string(k)=>v for (k,v) in pairs(p))
    d["time_vector"]=collect(0:p.Δt:p.T_move)
    d["sigma_eta_list"]=[0.0,0.019,0.04]
    d["random_seed"]=1; d["N_trials_MC"]=64
    for key in FIG2_TRAJECTORY_FIELDS; d[key]=Array{Float64,3}[]; end
    for key in ["K_TOD_all","L_TOD_all","K_LCO_all","L_LCO_all"]; d[key]=Array{Float64,3}[]; end
    for method in ["TOD","LCO"]
        for key in ["analytical_cost_","mean_MC_cost_","sem_MC_cost_"]; d[key*method]=Float64[]; end
        d["trial_costs_"*method]=Vector{Float64}[]
    end
    for σ_η in d["sigma_eta_list"]
        solution=solve_reaching(p,σ_η;N_iter_Todorov=3,N_iter_tot=3)
        for method in ["TOD","LCO"]
            K=getproperty(solution,Symbol("K_"*method));L=getproperty(solution,Symbol("L_"*method))
            mc=simulate_reaching(p,L,K,σ_η,64;random_seed=1)
            push!(d["K_"*method*"_all"],K);push!(d["L_"*method*"_all"],L)
            for state in ["x","z","u"]
                push!(d[state*"_store_"*method],permutedims(getproperty(mc,Symbol(state*"_traces")),(1,3,2)))
            end
            push!(d["analytical_cost_"*method],reaching_cost(p,L,K,σ_η))
            push!(d["mean_MC_cost_"*method],mc.empirical_cost)
            push!(d["sem_MC_cost_"*method],mc.empirical_cost_SE)
            push!(d["trial_costs_"*method],mc.cost_per_trial)
        end
    end
    d["validation_only"]=true;d["generated_by_current_core"]=true
    d["model"]=Dict(string(k)=>v for (k,v) in pairs(p))
    return d
end

function calibrated_reference()
    # Main Figure 2 is the single source for the model and internal-noise policies.
    data_dir=joinpath(@__DIR__,"..","..","..","data","published","Main","Figure2")
    reference=load_fig2(data_dir)
    p=reaching_problem()
    for (key,value) in pairs(p)
        isequal(reference[string(key)],value) || error("Figure 2 model mismatch: $key")
    end
    reference["sigma_eta_list"]==round.(range(0.0,stop=0.04,length=20),digits=3) ||
        error("The calibrated internal-noise reference must use the Figure 2 grid")
    reference["model"]=Dict(string(k)=>v for (k,v) in pairs(p))
    return reference
end

include("sensory_overestimation.jl")
include("sensory_underestimation.jl")

function simulate_identifiability(output_dir,experiments;validation_only=false)
    reference=validation_only ? validation_reference() : calibrated_reference()
    modules=Dict("overestimation"=>Overestimation,
        "underestimation"=>Underestimation)
    for experiment in experiments
        modules[experiment].run_experiment(reference,output_dir;validation_only)
    end
end
end
