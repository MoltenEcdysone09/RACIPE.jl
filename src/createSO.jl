using PackageCompiler

create_sysimage(
    [:CSV, :DataFrames, :Distributions, :ProgressMeter, :StatsBase, :OrdinaryDiffEq, :DiffEqCallbacks],
    sysimage_path="RACIPE.so",
    precompile_execution_file="RACIPE.jl" # the new line
)
