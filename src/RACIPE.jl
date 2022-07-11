module RACIPE

# Write your package code here.

include("read_topo.jl")
include("simul.jl")

export createRxnNet, rxnParamsList, genPrsFile, genParams, runRACIPE

end
