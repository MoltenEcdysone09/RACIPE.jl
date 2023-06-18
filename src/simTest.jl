println("IncRACIPE")
@time include("RACIPE.jl")
println("UseRACIPE")
@time using .RACIPE

println(Threads.nthreads())

println("CreateRxn")
@time createRxnNet("EMT.topo")

println("genPrs")
@time genPrsFile("EMT.topo"; numParas=1000)
@time genPrsFile("EMT.topo"; numParas=1000)

println("genParams")
@time genParams("EMT.prs"; numParas=1000)
@time genParams("EMT.prs"; numParas=1000)

println("IncRxnNet")
@time rxn = include("EMT.jl")

@time soldf = runRACIPE(rxn, "EMT_parameters.dat"; paramSets=1:1000, num_ini=100)
#println(first(soldf))
