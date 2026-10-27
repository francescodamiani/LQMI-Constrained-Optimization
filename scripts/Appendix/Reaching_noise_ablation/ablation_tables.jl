# Complete published CSV; no CSV.jl dependency or duplicate numerical exports.
using DelimitedFiles

const REACHING_ABLATION_CSV="calibrated_fixed_gain_noise_ablation_subset.csv"
const REACHING_ABLATION_RUN="calibrated_reaching_ablation.jld2"
const REACHING_ABLATION_CONDITIONS=["full_noise","no_process_noise",
    "no_additive_sensory_noise","no_all_sensory_noise","no_motor_multiplicative_noise",
    "no_sensory_multiplicative_noise","no_internal_noise"]
const REACHING_ABLATION_LEVELS=[0.0,0.011,0.019,0.029,0.040]

function load_reaching_ablation(data_dir)
    generated=joinpath(data_dir,REACHING_ABLATION_RUN)
    isfile(generated) && return JLD2.load(generated)
    original=joinpath(data_dir,REACHING_ABLATION_CSV)
    isfile(original) || error("Saved reaching ablation missing in $data_dir; use --simulate explicitly")
    values,header=readdlm(original,',',Any,'\n';header=true)
    names=Symbol.(vec(header))
    expected=[:sigma_eta,:condition,:terminal_elbow_SD_pct,:mean_angle_SD,:mean_velocity_SD,
        :terminal_elbow_SD_reduction_pct,:mean_angle_SD_reduction_pct,:mean_velocity_SD_reduction_pct]
    names==expected || error("Unexpected reaching ablation CSV columns")
    df=DataFrame([name==:condition ? String.(values[:,i]) : Float64.(values[:,i])
        for (i,name) in enumerate(names)],names)
    nrow(df)==35 || error("The reaching ablation must contain all 35 rows")
    unique(df.sigma_eta)==REACHING_ABLATION_LEVELS || error("Unexpected internal-noise grid")
    for σ_η in REACHING_ABLATION_LEVELS
        df.condition[df.sigma_eta.==σ_η]==REACHING_ABLATION_CONDITIONS || error("Incomplete ablation conditions")
    end
    return Dict("df_results"=>df,"generated_by_current_core"=>false)
end

function reaching_ablation_table(df,metric)
    return DataFrame([begin
        group=df[df.sigma_eta.==σ_η,:]
        values=[only(group[group.condition.==condition,metric]) for condition in REACHING_ABLATION_CONDITIONS]
        (;sigma_eta=σ_η,(Symbol(condition)=>value for (condition,value) in zip(REACHING_ABLATION_CONDITIONS,values))...)
    end for σ_η in unique(df.sigma_eta)])
end

function print_reaching_ablation(data;io=stdout)
    df=data["df_results"]
    println(io,"Calibrated Figure 2 reaching — fixed full-model LCO gains; ",
        get(data,"generated_by_current_core",false) ? "analytical recomputation" : "published results")
    println(io,"No reoptimization. Combined sensory removal zeros additive omega and multiplicative rho together.")
    columns=["eta"=>:sigma_eta,"Full"=>:full_noise,"No process"=>:no_process_noise,
        "No additive sensory"=>:no_additive_sensory_noise,"No combined sensory"=>:no_all_sensory_noise,
        "No motor multiplicative"=>:no_motor_multiplicative_noise,
        "No sensory multiplicative"=>:no_sensory_multiplicative_noise,"No internal"=>:no_internal_noise]
    for (label,metric) in [
        ("Terminal elbow SD (% of 20-degree reach)",:terminal_elbow_SD_pct),
        ("Terminal elbow SD reduction (%)",:terminal_elbow_SD_reduction_pct),
        ("Mean angular SD (rad; mean across 3 joints and times 2:41)",:mean_angle_SD),
        ("Mean angular SD reduction (%)",:mean_angle_SD_reduction_pct),
        ("Mean angular-velocity SD (rad/s; mean across 3 joints and times 2:41)",:mean_velocity_SD),
        ("Mean angular-velocity SD reduction (%)",:mean_velocity_SD_reduction_pct)]
        println(io,"\n",label)
        print_benchmark_table(io,reaching_ablation_table(df,metric),columns)
    end
    println(io,"\nReductions are 100*(1-SD_ablated/SD_full) with the same frozen gains. They are not additive variance fractions.")
    println(io,"Analytical terminal SD is 8.63464–10.54364%; the Figure 2 MC sample SD endpoints are 8.55236–10.45508%.")
    println(io,"Reported reduction ranges exclude eta=0; its process reduction is 83.7861% and its internal reduction is zero.")
    println(io,"Tables display 12 significant digits; saved CSV/JLD2 retains full available precision.")
end
