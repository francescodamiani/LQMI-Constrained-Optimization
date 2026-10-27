# Repository paths and output protection shared by figure scripts.
using Pkg, Dates
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)

function fig1_options(args)
    simulate = "--simulate" in args
    smoke = "--smoke" in args
    force = "--force" in args
    tag = ""
    data = ""
    panels = collect("bcdefg")
    for arg in args
        if startswith(arg, "--tag=")
            tag = split(arg, "="; limit=2)[2]
        elseif startswith(arg, "--data=")
            data = split(arg, "="; limit=2)[2]
        elseif startswith(arg, "--panels=")
            panels = collect(replace(split(arg, "="; limit=2)[2], "," => ""))
        elseif !(arg in ["--simulate", "--smoke", "--force", "--help"])
            error("Unknown argument: $arg")
        end
    end
    smoke && !simulate && error("--smoke requires --simulate")
    simulate && !isempty(data) && error("--data and --simulate cannot be combined")
    isempty(panels) && error("No panels selected")
    all(p -> p in "bcdefg", panels) || error("Panels must be among b,c,d,e,f,g")
    length(unique(panels)) == length(panels) || error("Repeated panel")
    panels = [p for p in "bcdefg" if p in panels]
    isempty(tag) || occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$", tag) || error("Invalid run tag")
    isempty(tag) && (tag = (smoke ? "smoke_" : "run_") * Dates.format(now(), "yyyymmdd_HHMMSSsss"))
    data_dir = isempty(data) ? joinpath(PROJECT_ROOT, "data", "published", "Main", "Figure1") :
               (isabspath(data) ? data : joinpath(PROJECT_ROOT, data))
    simulate && (data_dir = joinpath(PROJECT_ROOT, "data", "generated", "Main", "Figure1", tag))
    figure_dir = joinpath(PROJECT_ROOT, "figures", "generated", "Main", "Figure1")
    return (; simulate, smoke, force, tag, data_dir=normpath(data_dir), figure_dir, panels)
end

# Reject symlinks in output ancestry, including the generated root itself.
# --force can replace named files only inside these generated directories.
function prepare_generated_directory(path, kind; force=false)
    base = joinpath(PROJECT_ROOT, kind, "generated")
    relative = relpath(path, base)
    (relative == ".." || startswith(relative, "../") || isabspath(relative)) && error("Output must be under $base")
    current = PROJECT_ROOT
    for part in splitpath(relpath(path, PROJECT_ROOT))
        current = joinpath(current, part)
        islink(current) && error("Symlink output path is forbidden: $current")
    end
    ispath(path) && !isdir(path) && error("Output path is not a directory: $path")
    isdir(path) && !isempty(readdir(path)) && !force && error("Output exists: $path; choose another --tag or use --force")
    if isdir(path)
        for (root, dirs, files) in walkdir(path), entry in vcat(dirs, files)
            target = joinpath(root, entry)
            islink(target) && error("Symlink in output directory is forbidden: $target")
            isfile(target) && stat(target).nlink > 1 && error("Hard-linked output is forbidden: $target")
        end
    end
    mkpath(path)
    return path
end
