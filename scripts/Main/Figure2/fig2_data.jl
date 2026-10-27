# Complete calibrated archive: metadata/gains and one trajectory file per noise level.
using JLD2, CodecZlib, Printf

const FIG2_TRAJECTORY_FIELDS = ["$(state)_store_$(method)" for state in ["x","z","u"] for method in ["TOD","LCO"]]
fig2_trajectory_file(σ_η) = "trajectories_sigma_eta_" * replace(@sprintf("%.3f", σ_η), "." => "p") * ".jld2"

function save_fig2(data_dir, results)
    mkpath(data_dir)
    JLD2.jldopen(joinpath(data_dir, "sweep.jld2"), "w"; compress=true) do file
        for key in sort(collect(keys(results)))
            key in FIG2_TRAJECTORY_FIELDS || (file[key] = results[key])
        end
    end
    for (idx, σ_η) in enumerate(results["sigma_eta_list"])
        JLD2.jldopen(joinpath(data_dir, fig2_trajectory_file(σ_η)), "w"; compress=true) do file
            for key in FIG2_TRAJECTORY_FIELDS
                file[key] = results[key][idx]
            end
        end
    end
    return results
end

function load_fig2(data_dir)
    metadata = joinpath(data_dir, "sweep.jld2")
    isfile(metadata) || error("Missing Figure 2 data: $metadata. Use --simulate to recompute explicitly.")
    results = JLD2.load(metadata)
    files = [joinpath(data_dir, fig2_trajectory_file(σ)) for σ in results["sigma_eta_list"]]
    all(isfile, files) || error("Incomplete Figure 2 trajectory archive in $data_dir")
    for key in FIG2_TRAJECTORY_FIELDS
        results[key] = Array{Float64,3}[]
    end
    for file in files
        trajectories = JLD2.load(file)
        for key in FIG2_TRAJECTORY_FIELDS
            push!(results[key], trajectories[key])
        end
    end
    return results
end
