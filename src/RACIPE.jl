module RACIPE

# Write your package code here.
export createRxnNet, genParams, runRACIPE

include("read_topo.jl")
include("simul.jl")

end
