# Historical appendix Figure 6: fixed LCO gains, matched-variance noise shapes.
using Random, Statistics, Printf

function simulate_nongaussian_arm(; problem_data=arm_problem(),
        σ_η_internal_noise=0.1, N_trials_MC=100_000,
        N_iter_Todorov=500,N_iter_tot=500,verbose=true,batch_size=500)
    p = problem_data
    solution = solve_arm(p,σ_η_internal_noise;N_iter_Todorov,N_iter_tot,verbose)
    noise_specs = [(:gaussian,"Gaussian"),(:studentt,"Student-t (nu=5)"),(:beta,"Beta(2,5)")]
    trial_totals_dict = Dict{String,Vector{Float64}}()
    mean_cost_dict = Dict{String,Float64}()
    sem_cost_dict = Dict{String,Float64}()
    base_seed = 1234
    for (j,(noise_type,noise_label)) in enumerate(noise_specs)
        println("Running MC with $noise_label noise...")
        # One persistent stream per distribution; batching retains all trial costs
        # and avoids allocating several GB of discarded trajectories at once.
        rng = ArmNoiseRNG(MersenneTwister(base_seed+10_000*j),noise_type)
        trial_totals = zeros(N_trials_MC)
        for first_trial in 1:batch_size:N_trials_MC
            last_trial = min(first_trial+batch_size-1,N_trials_MC)
            mc = simulate_arm(p,solution.L_Full_Solver,solution.K_Full_Solver,
                σ_η_internal_noise,last_trial-first_trial+1;rng)
            trial_totals[first_trial:last_trial] .= mc.cost_per_trial
        end
        trial_totals_dict[noise_label] = trial_totals
        mean_cost_dict[noise_label] = mean(trial_totals)
        sem_cost_dict[noise_label] = std(trial_totals)/sqrt(N_trials_MC)
        println(@sprintf("  mean total cost = %.6f | SEM = %.6f",mean_cost_dict[noise_label],sem_cost_dict[noise_label]))
    end
    analytical_cost = arm_analytical_cost(p,solution.L_Full_Solver,solution.K_Full_Solver,σ_η_internal_noise)
    settings = (; σ_η_internal_noise,N_trials_MC,N_iter_Todorov,N_iter_tot,base_seed,
        trajectory_indexing="dimension,time,trial")
    return (; problem_data,settings,solution,trial_totals_dict,mean_cost_dict,sem_cost_dict,analytical_cost)
end
