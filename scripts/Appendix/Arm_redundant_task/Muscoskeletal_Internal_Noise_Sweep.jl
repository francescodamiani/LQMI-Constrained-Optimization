# Internal-noise sweep: TOD vs LCO for the redundant stabilization task.
using MultivariateStats, StatsBase, Statistics, LinearAlgebra, Random, Printf

function simulate_internal_noise_sweep(; problem_data=arm_problem(),
        sigma_eta_list=round.(range(0.0,stop=0.4,length=20),digits=2),
        N_trials_MC=500, N_iter_Todorov=500, N_iter_tot=500, verbose=true)
    p = problem_data
    Nnoise = length(sigma_eta_list)
    K_TOD_all = Vector{Array{Float64,3}}(undef,Nnoise)
    L_TOD_all = similar(K_TOD_all)
    K_Full_Solver_all = similar(K_TOD_all)
    L_Full_Solver_all = similar(K_TOD_all)
    best_cost_list_all = Vector{Vector{Float64}}(undef,Nnoise)
    x_store = similar(K_TOD_all)
    z_store = similar(K_TOD_all)
    u_store = similar(K_TOD_all)
    x_store_TOD = similar(K_TOD_all)
    z_store_TOD = similar(K_TOD_all)
    u_store_TOD = similar(K_TOD_all)
    pc_dirs_list = Vector{Matrix{Float64}}(undef,Nnoise)
    pc_dirs_list_TOD = similar(pc_dirs_list)
    pc_var_list = zeros(3,Nnoise)
    pc_var_list_TOD = zeros(3,Nnoise)
    mean_cost_MM, sem_cost_MM = zeros(Nnoise), zeros(Nnoise)
    mean_cost_TOD, sem_cost_TOD = zeros(Nnoise), zeros(Nnoise)
    analytical_cost_LCO, analytical_cost_TOD = zeros(Nnoise), zeros(Nnoise)
    cost_per_trial_LCO = Vector{Vector{Float64}}(undef,Nnoise)
    cost_per_trial_TOD = similar(cost_per_trial_LCO)
    solver_diagnostics = Vector{Any}(undef,Nnoise)

    for idx in eachindex(sigma_eta_list)
        σ_η_internal_noise = sigma_eta_list[idx]
        println("\n--- Internal noise level $σ_η_internal_noise ---")
        solution = solve_arm(p,σ_η_internal_noise; N_iter_Todorov,N_iter_tot,verbose)
        K_TOD_all[idx], L_TOD_all[idx] = solution.K_TOD, solution.L_TOD
        K_Full_Solver_all[idx], L_Full_Solver_all[idx] = solution.K_Full_Solver, solution.L_Full_Solver
        best_cost_list_all[idx] = solution.best_cost_list
        solver_diagnostics[idx] = (; best_method=solution.best_method, best_ρ=solution.best_ρ,
            cost_init=solution.cost_init, best_cost=solution.cost_Full_Solver,
            time_TOD=solution.t_TOD,time_LCO=solution.t_LCO)

        # Same per-noise, per-method seed as the source; new simulator draw order differs.
        random_seed = 1234
        mc_TOD = simulate_arm(p,solution.L_TOD,solution.K_TOD,σ_η_internal_noise,N_trials_MC;
            rng=MersenneTwister(random_seed))
        mc = simulate_arm(p,solution.L_Full_Solver,solution.K_Full_Solver,σ_η_internal_noise,N_trials_MC;
            rng=MersenneTwister(random_seed))
        # Native core indexing: dimension,time,trial. Controls come directly from MC.
        x_store[idx], z_store[idx], u_store[idx] = mc.x_traces, mc.z_traces, mc.u_traces
        x_store_TOD[idx], z_store_TOD[idx], u_store_TOD[idx] = mc_TOD.x_traces, mc_TOD.z_traces, mc_TOD.u_traces
        mean_cost_MM[idx], sem_cost_MM[idx] = mc.empirical_cost, mc.empirical_cost_SE
        mean_cost_TOD[idx], sem_cost_TOD[idx] = mc_TOD.empirical_cost, mc_TOD.empirical_cost_SE
        cost_per_trial_LCO[idx], cost_per_trial_TOD[idx] = mc.cost_per_trial, mc_TOD.cost_per_trial
        analytical_cost_LCO[idx] = arm_analytical_cost(p,solution.L_Full_Solver,solution.K_Full_Solver,σ_η_internal_noise)
        analytical_cost_TOD[idx] = arm_analytical_cost(p,solution.L_TOD,solution.K_TOD,σ_η_internal_noise)

        # Pooled-control PCA and its first three fractions.
        for (u_vec, pc_dirs, pc_var) in [(mc.u_traces,pc_dirs_list,pc_var_list),
                                        (mc_TOD.u_traces,pc_dirs_list_TOD,pc_var_list_TOD)]
            Umat = reshape(u_vec,p.dimension_of_control,:)
            pca_model = fit(PCA,Umat;maxoutdim=p.dimension_of_control)
            eigvals = principalvars(pca_model)
            pc_dirs[idx] = projection(pca_model)
            pc_var[:,idx] .= (eigvals ./ sum(eigvals))[1:3]
        end
        println(@sprintf("TOD cost %.6g ± %.3g; LCO %.6g ± %.3g",mean_cost_TOD[idx],
            sem_cost_TOD[idx],mean_cost_MM[idx],sem_cost_MM[idx]))
    end
    settings = (; N_trials_MC,N_iter_Todorov,N_iter_tot,random_seed=1234,
        trajectory_indexing="dimension,time,trial")
    return (; problem_data,sigma_eta_list,settings,K_TOD_all,L_TOD_all,
        K_Full_Solver_all,L_Full_Solver_all,best_cost_list_all,
        x_store,z_store,u_store,x_store_TOD,z_store_TOD,u_store_TOD,
        pc_dirs_list,pc_dirs_list_TOD,pc_var_list,pc_var_list_TOD,
        mean_cost_MM,sem_cost_MM,mean_cost_TOD,sem_cost_TOD,
        analytical_cost_LCO,analytical_cost_TOD,cost_per_trial_LCO,cost_per_trial_TOD,
        solver_diagnostics)
end
