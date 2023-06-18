##
##using DifferentialEquations
#using CSV
#using DataFrames
#using ModelingToolkit
#using OrdinaryDiffEq
#using Distributions
#using ProgressMeter
#using StatsBase
##using DataStructures: SortedDict, Reverse
##using Combinatorics

#=
Function: round_log!
Class: side effect
Input: Float64 type number
Output: Input rounded to 6 digits. Can be changed to give log2 values (Testing)
Description: Takes in a number and rounds it to 6 decimal places. Is also capable of giving the log2 values of the steady state values (Needs to be tested wrt indetification of steady states).
=#
##
function round_log!(x::Float64)
    #if x < 0 && x > -1
    #    x = 1e-9
    #end
    #x = round(log2(x), digits = 2)
    x = round(x, digits=6)
    return x
end

#=
Function: prob_func
Class: side effect
Input: ODE Probelm  made from the the reaction network (prob)
Output: Ensemble poblem with different initial conditions.
Description: Takes in an already created ODE probelm to create an ensemble problem which is simulated over mutiple initial conditions.
=#
##
function prob_func(prob,i,repeat)
    remake(prob,u0=rand(Uniform(0,100),length(prob.u0)))
end

#=
Function: identifyStates!
Class: proper
Input: State count dictionary with the states (as list) as key and the corresponding counts as the values.
Output: Removes and merges duplicate steady states.
Description: Identifies all the states which are within euclidean distance threshold of 1 with some othre state and merges them. Uses proportionmap output as its input.
=#
##
function identifyStates!(st_counts::Dict{Any, Float64})
    st_key = collect(keys(st_counts))::Vector{Any}
    stk = [st_key[1]]::Vector{Vector{Float64}}
    for st in st_key
        dupli = false
        for r in stk
            if StatsBase.L2dist(st, r) < 1
                dupli = true
                if st != st_key[1]
                    st_counts[r] = st_counts[r] +  st_counts[st]
                    delete!(st_counts, st)
                end
                break
            end
        end
        if dupli == false
            push!(stk, st)
        end
    end
    return st_counts
end

#=
Function: roundStates!
Class: proper
Input: Simulation results (simu_results) of the ensemble simualtions over mutiple initial conditions.
Output: A dict with steady state (list) as the key and the corresponding proportion as the value.
Description: Takes in the raw simualtion results as input, rounds the steady states, and then passes it to proprtionmap function, it then gives the proprtion of the different steady states. This output will be taken in by identifyStates! to merge some duplicates still remaining to give final steady state frequencies.
=#
##
function roundStates!(simu_results, num_ini::Int64)
    sols = []
    for i in 1:num_ini
        sol = round_log!.(simu_results[i].u[end])
        append!(sols, [sol])
    end
    return proportionmap(sols)
end

#=
Function: genParamMatrix
Class: proper
Input: Paramter file name.
Output: A matrix of the paramters.
Description: Takes in  the paramter files, parses it into a dataframe, and then converts the dataframe into a matrix. Matrix because it was observed that it lead to a performance boost.
=#
##
function genParamMatrix(param_file::String)
    param_df = CSV.read(param_file, DataFrame)
    return Matrix(param_df)::Matrix{Float64}
end

#=
Function: genParamCols
Class: proper
Input: Paramter file name.
Output: Paramter names.
Description: Takes in  the paramter files, parses it into a dataframe, and then gives the columns as a list.
=#
##
function genParamCols(param_file::String)
    param_df = CSV.read(param_file, DataFrame)
    param_names = names(param_df)
    #node_names = [i[6:end] for i in param_names if findfirst("Prod_", i) == 1:5]
    node_names = [i[6:end] for i in param_names if occursin("Prod_", i)]
    print(node_names)
    return param_names, node_names
end

#=
Function: simulateEnsemble!
Class: proper
Input: paramter set (p), number of initial conditions (num_ini)
Output: Steady states and their frequencies for that parameter set.
Description: Takes in paramter set and the number of initial conditions as the input, creates an ODE Problem, which is then converted into a EnsembleProblem. The EnsembleProblem is then simulated over the initial condiations (parallel simulations), the output is then passed to indentify steady states function, which then gives the final output.
=#
##
function simulateEnsemble!(RxnNet, p, num_ini::Int64, num_nodes::Int64)
    u0 = rand(Uniform(0,100), num_nodes)::Vector{Float64}
    tspan = (0.0, 100.0)::Tuple{Float64, Float64}
    prob = ODEProblem{true}(RxnNet, u0, tspan, p)
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
    sim = solve(ensemble_prob, Tsit5(),save_everystep=false,save_start=false,callback=TerminateSteadyState(1e-3, 1e-5),EnsembleThreads(),trajectories=num_ini)
    return identifyStates!(roundStates!(sim, num_ini))
end

#=
Function: runRACIPE [Needs Improvement]
Class: proper
Input: The paramter file; para_range - A vector or a range (x:y) of the parameters to be simulated , number of initial conditions
Output: A dataframe with the paramters and thier corresponding steady states
Description: Parses the paramter file into a matrix and then converts it into a matrix. Loops through the matrix rows to simulate the network over multiple initial conditions. Finally, the steady states identified are then compiles into a dataframe along with thier paramter sets (similar to solution.dat).
=#
##
function runRACIPE(rxnNet, param_file::String; paramSets=1:1, num_ini::Int64=1000)
    #param_df = readParameters(param_file)
    #param_df = CSV.read(param_file, DataFrame)
    #node_names = [i[2:end] for i in param_names if findfirst("g", i) == 1:1]
    # Generates a paramter matrix
    param_df = genParamMatrix(param_file)
    # A tuple with paramters names and node names
    names = genParamCols(param_file)
    # paramter names list
    param_names = names[1]
    # node names list
    node_names = names[2]
    # get number of nodes in a reaction network
    num_nodes = length(node_names)
    # Initialise an empty vector to store all the solutions in the format
    # "ParamNo","ParameterNames"...,"Steady State Values"...,"Relabtive Stability"
    solMatrix = []::Vector{Any}
    # Loop througha all the parameters to solve them and get thier solutions
    @showprogress for f in paramSets
        simu_results = simulateEnsemble!(rxnNet, param_df[f,:], num_ini, num_nodes)
        for (key, value) in simu_results
            push!(solMatrix, vcat(f, param_df[f,:], key, round(value, digits=4)))
        end
    end
    # Convert the vector of vector into a matrix (' added at last to transpose)
    solMatrix = reduce(hcat,solMatrix)'
    # Create a column names list
    solCols = vcat("ParamNo",param_names,node_names,"RelStab")
    # Return the dataframe of solutions
    return DataFrame(solMatrix, solCols)
end

export runRACIPE

