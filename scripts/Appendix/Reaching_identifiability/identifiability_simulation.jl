# Shared adapters to the validated library; no duplicated stochastic dynamics.
module IdentifiabilitySimulation
using LinearAlgebra, Random, Statistics, JLD2
include(joinpath(@__DIR__,"..","..","..","src","MyFunctionsOptimalRegularizedSolutions.jl"))
include(joinpath(@__DIR__,"..","..","Main","Figure2","fig2_data.jl"))
include(joinpath(@__DIR__,"..","..","Main","Figure2","reaching_task.jl"))

function sample_identifiability(seed, N, n, nz, p, T, mx, mz, Σx, Σz,
    A, B, C, H, D, K, M, L, Ωξ, Ωω, Ωη, Q, R)
    # Historical scripts expect (dimension, trial, time); the core returns
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
    # Same dynamics/cost construction, with the historical .01 external noises.
    # Only horizon, trials, iterations and sweep sizes are reduced for validation.
    p0 = reaching_problem(T_move=0.07)
    p = merge(p0,(;σ_ω_add_sensory_noise=0.01,σ_ϵ_mult_control_noise=0.01,
        σ_ρ_mult_sensory_noise=0.01,C=0.01*p0.B,D=0.01*p0.H,
        Ω_ω=Matrix(Diagonal(fill(0.01^2,3)))))
    d = Dict{String,Any}(string(k)=>v for (k,v) in pairs(p))
    d["time_vector"]=collect(0:p.Δt:p.T_move)
    d["sigma_eta_list"]=[0.0,0.01,0.03]
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
    return d
end

include("sensory_overestimation.jl")
include("matched_noise_sources.jl")
include("sensory_underestimation.jl")

function simulate_identifiability(output_dir,published_dir,experiments;validation_only=false)
    reference = validation_only ? validation_reference() :
        (any(in(["overestimation","noise_sources"]),experiments) ?
            load_fig2(joinpath(published_dir,"reference")) :
            JLD2.load(joinpath(published_dir,"reference","sweep.jld2")))
    if validation_only
        save_fig2(joinpath(output_dir,"reference"),reference)
    else
        # Record the model/reference gains without duplicating archived trajectories.
        mkpath(joinpath(output_dir,"reference"))
        JLD2.jldopen(joinpath(output_dir,"reference","sweep.jld2"),"w";compress=true) do file
            for (key,value) in reference
                key in FIG2_TRAJECTORY_FIELDS || (file[key]=value)
            end
        end
    end
    modules=Dict("overestimation"=>Overestimation,"noise_sources"=>MatchedNoiseSources,
        "underestimation"=>Underestimation)
    for experiment in experiments
        modules[experiment].run_experiment(reference,output_dir;validation_only)
    end
end
end
