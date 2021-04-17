using LightGraphs, SimpleWeightedGraphs
using GraphPlot
using BlockArrays 
using LinearAlgebra
using PinningControl 

n = 5 
g1 = cycle_graph(n) 
g2 = cycle_graph(n) 
g3 = blockdiag(g1, g2) |> SimpleWeightedDiGraph 

add_edge!(g3, 5, 10, -1.)
add_edge!(g3, 10, 5, -1)
add_edge!(g3, 5, 6, 1.)
add_edge!(g3, 6, 5, 1.)

gplot(g3, nodelabel=1:nv(g3)) |> display
Φ = -collect(laplacian_matrix(g3))
Φ = BlockMatrix(Φ, [n, n], [n, n])
display(Φ)
issymmetric(Φ) 
is_strongly_connected(SimpleDiGraph(getblock(Φ, 1, 1)))
is_strongly_connected(SimpleDiGraph(getblock(Φ, 2, 2)))
