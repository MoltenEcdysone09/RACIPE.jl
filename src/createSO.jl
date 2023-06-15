using PackageCompiler

create_sysimage(
    #[:CSV, :DataFrames, :DelimitedFiles, :Distributions, :OrderedCollections, :Statistics, :DifferentialEquations, :ProgressMeter, :StatsBase, :ModelingToolkit, :OrdinaryDiffEq],
    [:CSV, :DataFrames, :Distributions, :OrderedCollections, :DifferentialEquations, :ProgressMeter, :StatsBase, :ModelingToolkit, :OrdinaryDiffEq],
    sysimage_path="RACIPE.so",
    precompile_execution_file="RACIPE.jl" # the new line
)
