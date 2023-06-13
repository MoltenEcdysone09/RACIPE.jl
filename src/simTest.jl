println("IncRACIPE")
@time include("RACIPE.jl")
println("UseRACIPE")
@time using .RACIPE

println(Threads.nthreads())

println("CreateRxn")
@time createRxnNet("EMT.topo")

println("genPrs")
@time genPrsFile("EMT.topo"; numParas=10000)
@time genPrsFile("EMT.topo"; numParas=10000)

println("genParams")
@time genParams("EMT.prs"; numParas=10000)
@time genParams("EMT.prs"; numParas=10000)

println("IncRxnNet")
@time rxn = include("EMT.jl")

#@time soldf = runRACIPE(rxn, "EMT_parameters.dat"; paramSets=1:10000, num_ini=1)
@time runRACIPE(rxn, "EMT_parameters.dat"; paramSets=1:2, num_ini=100)
@time soldf = runRACIPE(rxn, "EMT_parameters.dat"; paramSets=1:10000, num_ini=100)
#@time soldf = runRACIPE(rxn, "EMT_parameters.dat"; paramSets=1:10000, num_ini=10000)
#println(soldf)
