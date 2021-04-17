# This file includes an example of full synchronization in undirected cluster synchronization 

using PinningControl 
using Plots 
using LightGraphs    
using LinearAlgebra

# Settings 
node = Chua() 
θ = 10 
graph = random_regular_graph(100, 5) 
# graph = star_graph(10) 
Ξ = -collect(laplacian_matrix(graph))
# Ξ = [
#     -1 1 0 0 0 0; 
#     1 -3 1 1 0 0; 
#     0 1 -1 0 0 0; 
#     0 1 0 -3 1 1; 
#     0 0 0 1 -1 0; 
#     0 0 0 1 0 -1
# ]
n = size(Ξ, 1)
d = 3 
l = n ÷ 4
J = 1 : l
P = I(d)
Δ = 10 * I(d)


# Construct and solve network
net = undirectedfull(node, Ξ, P, Δ, J)
sol = netsolve(net, 0., 0.01, 1.)

# Plot solution 
t = sol.t 
x = sol.u 
plt1 = plot(getindex.(x, 1), getindex.(x, 2))
nr = 1 : 100
tr = t[nr] 
xr = x[nr] 
sr = getindex.(xr, 1)
plt2 = plot()
for k in d + 1 : d : n * d + d
    plot!(tr, sr - getindex.(xr, k))
end 

display(plot(plt1, plt2, layout=(2,1)))
