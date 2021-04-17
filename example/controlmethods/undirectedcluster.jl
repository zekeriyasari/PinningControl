# This file includes and example script for cluster synchronization in undirected networks 

using PinningControl 
using Plots 
using LightGraphs    
using LinearAlgebra
using GraphPlot

# Settings 
node = Chua()
n = 6
J = [n ÷ 2, n] 
Φ = initmat(J)
graph = SimpleGraph(Φ)
gplot(graph, nodelabel=1:nv(graph))
# Φ = [
#     -1. 1 0 0 0 0;
#     1 -2 1 0 -1 1;
#     0 1 -1 0 1 -1;
#     0 0 0 -1 1 0;
#     0 -1 1 1 -2 1;
#     0 1 -1 0 1 -1
#     ]
# n = size(Φ, 1) 
d = 3 
# J = [n ÷ 2, n]
l = length(J) 
P = I(d) 
Δ = 10 * I(d)
γ = 200 * ones(l)

# Construct the network 
net = undirectedcluster(node, Φ, γ, Δ, P, J, rand(n * d))
sol = netsolve(net, 0., 0.01, 10.)

# Plot solution 
t = sol.t 
x = sol.u 

plt1 = plot(getindex.(x, 1), getindex.(x, 2))

nr = 1 : 100
tr = t[nr] 
xr = x[nr] 
ks = PinningControl.clusters(J)

clsidx = 1
sr1 = getindex.(xr, (clsidx - 1) * d + 1)
plt2 = plot()
for k in ks[clsidx] .+ l
    xi =  getindex.(xr, (k - 1) * d + 1) 
    plot!(tr, sr1 - xi,  label="$k")
end 

clsidx = 2
sr2 = getindex.(xr, (clsidx - 1) * d + 1)
plt3 = plot()
for k in ks[clsidx] .+ l
    xi =  getindex.(xr, (k - 1) * d + 1) 
    plot!(tr, sr2 - xi,  label="$k")
end 

plt4 = plot(t, abs.(getindex.(x, 1) - getindex.(x, 4)))

display(plot(plt1, plt2, plt3, plt4))
