module PinningControl

using DifferentialEquations 
using LightGraphs, SimpleWeightedGraphs
using BlockArrays 
using LinearAlgebra

include("dynamics.jl") 
include("network.jl") 
include("pinningcontrol.jl") 
include("utils.jl") 

export 
Lorenz, Chua, Chen, 
undirectedfull, undirectedcluster, directedfull, directedcluster,
netprob, netsolve,
isirreducible, isdiffusive, isnondiffusive, initmat

end # module
