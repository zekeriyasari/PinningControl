# This file includes utility functions 

isirreducible(mat::AbstractMatrix) = is_strongly_connected(SimpleDiGraph(mat))


isdiffusive(mat::AbstractMatrix) = iszero(sum(mat, dims=2))


isnondiffusive(mat::AbstractMatrix) = mat[[idx for idx in 1 : length(mat) if idx ∉ diagind(mat)]] .≥ 0


# Helper functions for construction of initial connection matrices 

function initmat(J::AbstractVector) 
    l = length(J) 
    ks = length.(clusters(J))
    subgraphs = map(k -> cycle_graph(k), ks) 
    graph = blockdiag(subgraphs[1], subgraphs[2])
    for n in 3 : length(subgraphs)
        graph = blockdiag(graph, subgraphs[n])
    end 
    graph = SimpleWeightedDiGraph(graph) 
    for i in 1 : length(J) - 1 
        add_edge!(graph, J[i], J[i + 1], -1.) 
        add_edge!(graph, J[i + 1], J[i], -1.) 
        add_edge!(graph, J[i], J[i] + 1, 1.) 
        add_edge!(graph, J[i] + 1, J[i], 1.) 
    end 
    Φ = -collect(laplacian_matrix(graph))
    all(sum(Φ, dims=2) .== 0) || @warn "Row sum of the matrix is not all zero"
    issymmetric(Φ) || @warn "The matrix is not symmetric"
    Φ
end 

# function initmat(J::AbstractVector, n::Int, topology=complete_graph)
#     l = length(J) 
#     Gs = clusters(J)
#     k = length.(Gs)
#     Φ = BlockArray(zeros(n, n), k, k) 
#     for i in 1 : l 
#         setblock!(Φ, ondiagmat(topology, k[i]), i, i)
#     end 
#     for i in 1 : l - 1, j in i + 1 : l 
#         bmat = offdiagmat(k[i], k[j])
#         setblock!(Φ, bmat, i, j)
#         setblock!(Φ, bmat', j, i)
#     end 
#     # all(sum(Φ, dims=2) .== 0.) || error("Zero row sum property is not satisfied\nΦ:$(Φ)")
#     Φ
# end

function ondiagmat(topology=complete_graph, args...; kwargs...) 
    graph = topology(args...; kwargs...)
    is_strongly_connected(SimpleDiGraph(graph)) ||  error("Could not costruct irreducible matrix")
    mat = -collect(laplacian_matrix(graph))
    mat 
end 

function offdiagmat(nrows, ncols) 
    mat = rand(0.: 1., nrows, ncols) 
    idx = diagind(mat) 
    mat[idx] .= 0 
    mat[idx] -= sum(mat, dims=2)
    (mat + mat') / 2
end 
