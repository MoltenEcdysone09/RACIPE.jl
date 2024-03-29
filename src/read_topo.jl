#=
Created by: Pradyumna
Created on: 18/08/2020
Brief description: Functions to create a reaction network and the corresponding paramters from a topo file.
=#
#using DelimitedFiles
#using DataFrames
#using Distributions
#using OrderedCollections
#using Statistics
#using DifferentialEquations
#using CSV


#=
Function: node_add
Class: side_effect (i.e., no output)
Input: topo_row -> a row in the topo file
Output: node_dict (this is not in the function environment. should it be supplied as an argument?)
Description: Adds pr_node to node_dict such that the Source and type of the Target node in pr_node are listed under its key in node_dict.
=#
function node_add(topo_row::SubArray, node_dict::Dict)
    if haskey(node_dict, topo_row[2]) == false
        #print(node_dict[pr_node])
        node_dict[topo_row[2]] = [topo_row[1], topo_row[3]]
    else
        push!(node_dict[topo_row[2]], topo_row[1], topo_row[3])
    end
end


#=
Function: net_dict
Class: proper (i.e., outputs a variable)
Input: topo_dat -> 2-d array of type Any with 3 columns, each row lists an edge in the network in the form of "Source Target Type".
Output: node_dict
Description: Takes the array generated from topo file and returns a dictionary that, for each node in the network (key) lists the incoming edges (value)
=#
function net_dict_gen(topo_dat::Array{Any})
    node_dict = Dict()
    for rw in eachrow(topo_dat)
        node_add(rw, node_dict)
    end
    #= Addtion of nodes which have only out-going edges =#
    for rw in eachrow(topo_dat)
        #println(rw)
        if haskey(node_dict, rw[1]) == false
            #println(rw[1])
            node_dict[rw[1]] = [rw[1], 0]
        end
    end
    return node_dict
end

#=
Function: stdy_st
Class: proper (i.e., outputs a variable)
Input: Array with [n0(value to be update), M0 (Median value for sampling), edge_type).
Output: Returns updated steady state value of n0
Description: Takes an array and then generates random parameters (with M0 for sampling for th) and then returns an updated steady state value of n0 .
=#
function stdy_st(n0::Float64, M0::Float64, edg_typ::Int)
    g = rand(Uniform(1, 100))
    k = rand(Uniform(0.1, 1))
    th = rand(Uniform(0.02 * M0, 1.98 * M0))
    n = round(rand(Uniform(1, 6)))
    f = rand(Uniform(1, 100))
    if edg_typ == 1
        return n0 * (f + (1 - f) * ((th^n) / (th^n + (g / k)^n))) / f
    else
        f = 1/f
        return n0 * (f + (1 - f) * ((th^n) / (th^n + (g / k)^n)))
    end
end

#=
Function: thr_gen
Class: proper (i.e., outputs a variable)
Input: initial parameters, number of Ativations, number of Inhibitions
Output: Median of thr_ary which gives the seed for sampling threshold
Description: Takes the number of initial paramters, number of activations and number of inhibitions; uses stdy_st to update the thr_ary from which the median value is returned, this value will be used for threshold samling.
=#
function thr_gen(init_param::Int, numA::Int, numI::Int)
    #Generate M0
    M0 = median(
                rand(Uniform(1, 100), init_param) ./ rand(Uniform(0.1, 1), init_param),
               )
    #Create an empty threshold array
    thr_ary = []
    #Start Loop to find out the final threshold
    for t = 1:init_param
        #Add g/k to the array, this will be inserted at thr_ary[t]
        push!(thr_ary, rand(Uniform(1, 100)) / rand(Uniform(0.01, 1)))
        #Update this value (thr_ary[t]) accroding to the nature of the edge
        if numA != 0
            for r = 1:numA
                thr_ary[t] = stdy_st(thr_ary[t], M0, 1)
            end
        end
        if numI != 0
            for r = 1:numI
                thr_ary[t] = stdy_st(thr_ary[t], M0, 2)
            end
        end
    end
    return median(thr_ary)
end

#=
Function: createRxnNet
Class: proper (i.e., outputs a variable) + side effect
Input: Name of the topology file (*.topo)
Output: List of parameters (param_list)  + A reaction network file representing the network through a systems of ODEs
Description: Takes the network topology as an input through the topo file, parses the file and create a .jl file of the same name havig the coupled odes representing the topology as the output file. Also reutrns the paramter list which is generted during this process, which will be used by param_gen function to generate paramters.
=#
function createRxnNet(topoFile::String)
    # Create the reaction network file
    rnfl = open(replace(topoFile, "topo" => "jl"), "w")
    #write Hsn and Hsp ration functions in reation network file
    #write shifted hill equation for negetive regulation
    write(rnfl, "#Defining custom Shifted Hill Equations\n")
    write(rnfl, "#Shifted Hill Equation for -ve regulation\n")
    write(
          rnfl,
          "Hsn(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))\n",
         )
    #write shifted hill equation for positive regualtion
    write(rnfl, "#Shifted Hill Equation for +ve regulation\n")
    write(
          rnfl,
          "Hsp(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))/Fld\n\n",
         )
    # Intialise the parameter lists
    param_list = String[]
    # Production paramters
    paramg_list = String[]
    # Degredation paramters
    paramk_list = String[]
    # List of the lines that need to be writen in the reaction file
    rxn_lines = String[]
    # Read the topo file as a delimited file
    #topo_dat = readdlm(topoFile)
    topo_dat = CSV.read(topoFile, DataFrame) 
    sort!(topo_dat, [:Target, :Source, :Type])
    topo_dat = Matrix(topo_dat)
    # Removing the empty line
    #topo_dat = topo_dat[1:end.!=1, :]
    # Get the dictionary of nodes : their edge and type of edge
    node_dict = net_dict_gen(topo_dat)
    #Start writing the reaction netowrk in a julia file.
    write(rnfl, "#Defining the reaction network with ModellingToolkit.jl macro\n")
    write(rnfl, "function rn!(du, u, p, t)\n")
    # Counter for number of the differntial equation
    to = 1
    # Loop to add the different reactions node wise in the ModellingToolkit.jl macro format
    for k in keys(node_dict)
        shill_parts = []
        push!(shill_parts, "Prod_$k")
        push!(paramg_list, "Prod_$k")
        push!(paramk_list, "Deg_$k")
        for n = 1:2:length(node_dict[k])
            cn = convert(String, node_dict[k][n])
            if node_dict[k][n+1] == 2
                push!(shill_parts, "Hsn(InhThr_$cn"*"_$k"*",$cn"*",InhFld_$cn"*"_$k"*",Hill_$cn"*"_$k"*")")
                push!(param_list, "InhThr_$cn"*"_$k", "InhFld_$cn"*"_$k", "Hill_$cn"*"_$k")
            elseif node_dict[k][n+1] == 1
                push!(shill_parts, "Hsp(ActThr_$cn"*"_$k"*",$cn"*",ActFld_$cn"*"_$k"*",Hill_$cn"*"_$k"*")")
                push!(param_list, "ActThr_$cn"*"_$k", "ActFld_$cn"*"_$k", "Hill_$cn"*"_$k")
            end
        end
        push!(rxn_lines, "\t" * "du[" * string(to) * "] = " * join(shill_parts, "*") * " - Deg_$k"*"*$k"*"\n")
        to = to + 1
    end
    #
    write(rnfl, "\t" * join(keys(node_dict), ", ") * " = u\n")
    #= Write all paramteres in the order, the same order needs to be followed in the
    paramteres file =#
    write(
          rnfl,
          "\t" *
          join(paramg_list, ", ") *
          " ," *
          join(paramk_list, ", ") *
          " ," *
          join(param_list, ", ") *
          " = p\n"
         )
    # Write all the lines in rxn_lines into the file
    for l in 1:length(rxn_lines)
        write(rnfl, rxn_lines[l])
    end
    # End the reaction function in the rxn file
    write(rnfl, "end\n")
    # Close the reaction file
    close(rnfl)
    # Return the list of paramters
    #return param_list = vcat(paramg_list, paramk_list, param_list)
end

#=
Function: rxnParamsList
Class: proper
Input: Name of the topology file (*.topo)
Output: Paramter Names List
Description: Generates paramters list corresponding to the topology file
=#
##
function rxnParamsList(topoFile::String)
    # Intialise the parameter lists
    param_list = String[]
    # Production paramters
    paramg_list = String[]
    # Degredation paramters
    paramk_list = String[]
    # List of the lines that need to be writen in the reaction file
    rxn_lines = String[]
    # Read the topo file as a delimited file
    #topo_dat = readdlm(topoFile)
    topo_dat = CSV.read(topoFile, DataFrame) 
    sort!(topo_dat, [:Target, :Source, :Type])
    topo_dat = Matrix(topo_dat)
    # Removing the empty line
    #topo_dat = topo_dat[1:end.!=1, :]
    # Get the dictionary of nodes : their edge and type of edge
    node_dict = net_dict_gen(topo_dat)
    # Counter for number of the differntial equation
    to = 1
    # Loop to add the different reactions node wise in the ModellingToolkit.jl macro format
    for k in keys(node_dict)
        shill_parts = []
        push!(shill_parts, "Prod_$k")
        push!(paramg_list, "Prod_$k")
        push!(paramk_list, "Deg_$k")
        for n = 1:2:length(node_dict[k])
            cn = convert(String, node_dict[k][n])
            if node_dict[k][n+1] == 2
                push!(param_list, "InhThr_$cn"*"_$k", "InhFld_$cn"*"_$k", "Hill_$cn"*"_$k")
            elseif node_dict[k][n+1] == 1
                push!(param_list, "ActThr_$cn"*"_$k", "ActFld_$cn"*"_$k", "Hill_$cn"*"_$k")
            end
        end
        to = to + 1
    end
    # Return the list of paramters
    return param_list = vcat(paramg_list, paramk_list, param_list)
end


#=
Function: genPrsFile
Class: side-effect
Input: Name of the topology file (*.topo), List of the paramters (same as that in the rxn file), number of parameters
Output: Paramter File
Description: Generates paramters file corresponding to the topology file which contains their Min and Max values.
=#
##
function genPrsFile(topoFile::String; numParas::Int64 = 1000)
    #Create a lsit of the parameters from the topo file
    param_list = rxnParamsList(topoFile)
    # Read the topo file as a delimited file
    #topo_dat = readdlm(topoFile)
    topo_dat = CSV.read(topoFile, DataFrame) 
    sort!(topo_dat, [:Target, :Source, :Type])
    topo_dat = Matrix(topo_dat)
    # Removing the empty line
    #topo_dat = topo_dat[1:end.!=1, :]
    # Get the dictionary of nodes : their edge and type of edge
    node_dict = net_dict_gen(topo_dat)
    #Create an empty dataframe for the parameters
    param_df = DataFrame( Parameter = String[], MinValue = Float64[], MaxValue = Float64[])
    # Populating the parameter dataframe with generated parameter values
    for p in param_list
        if occursin("Prod_", p)
            push!(param_df, [p, 1, 100])
        elseif occursin("Deg_", p)
            push!(param_df, [p, 0.1, 1])
        elseif occursin("InhThr_", p)
            push!(param_df, [p, 0, 0])
        elseif occursin("ActThr_", p)
            push!(param_df, [p, 0, 0])
        elseif occursin("Hill_", p)
            push!(param_df, [p, 1, 6])
        elseif occursin("InhFld_", p)
            push!(param_df, [p, 0.1, 1])
        elseif occursin("ActFld_", p)
            push!(param_df, [p, 1, 100])
        end
    end
    # Genrating threshold values and replacing them in the datafame
    for k in keys(node_dict)
        for n = 1:2:length(node_dict[k])
            con_nd = node_dict[node_dict[k][n]]
            numA = 0
            numI = 0
            for n = 2:2:length(con_nd)
                if con_nd[n] == 1
                    numA = numA + 1
                elseif con_nd[n] == 2
                    numI = numI + 1
                end
            end
            th = thr_gen(numParas, numA, numI)
            if node_dict[k][n+1] == 2
                param_df[param_df.Parameter .== "InhThr_"*node_dict[k][n]*"_"*k, :MinValue] .= 0.02 * th
                param_df[param_df.Parameter .== "InhThr_"*node_dict[k][n]*"_"*k, :MaxValue] .= 1.98 * th
            elseif node_dict[k][n+1] == 1
                param_df[param_df.Parameter .== "ActThr_"*node_dict[k][n]*"_"*k, :MinValue] .= 0.02 * th
                param_df[param_df.Parameter .== "ActThr_"*node_dict[k][n]*"_"*k, :MaxValue] .= 1.98 * th
            end
        end
    end
    # Write the dataframe in a csv file
    CSV.write(replace(topoFile, "topo" => "prs"), param_df; delim="\t")
end

#=
FunctionParams
Class: side effect
Input: Parameter file (*.prs); number of parameters
Output: Paramter Set Dataframe
Description: Generates paramters sets corresponding to the ranges given in the paramter file (*.prs).
=#
function genParams(prsFile::String; numParas::Int64 = 1000)
    # Read the topo file as a delimited file
    prsdf = CSV.read(prsFile, DataFrame)
    # Parameter column is the parameter list
    param_list = prsdf[!,:Parameter]
    #Create an empty dataframe for the parameters
    param_df = DataFrame()
    # Populating the parameter dataframe with generated parameter values
    for p in param_list
        minval = prsdf[prsdf.Parameter .== p, :MinValue][1]
        maxval = prsdf[prsdf.Parameter .== p, :MaxValue][1]
        if occursin("Hill_", p)
            param_df[!, p] = round.(rand(Uniform(minval, maxval), numParas))
        else
            param_df[!, p] = rand(Uniform(minval, maxval), numParas)
        end
    end
    # Write the dataframe in a csv file
    CSV.write(replace(prsFile, ".prs" => "_parameters.dat"), param_df)
end

export createRxnNet, genParams, rxnParamsList, genPrsFile

