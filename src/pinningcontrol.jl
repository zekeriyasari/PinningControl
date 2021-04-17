# This file includes 

function undirectedfull(node::AbstractNodeDynamics, 
                        Ξ::AbstractMatrix, 
                        P::AbstractMatrix, 
                        Δ::AbstractMatrix, 
                        J::AbstractVector, 
                        x0::AbstractVector=rand(size(Ξ,1)*size(P,1))*1e-3)
    l = length(J)
    d = size(P, 1)
    n = size(Ξ, 1)
    θ = maximum(eigvals((Δ + Δ'))) / 2
    Ξ22 = Ξ[l + 1 : end, l + 1 : end]
    ϵ = θ / abs(maximum(eigvals(Ξ22))) * 2  # The threshold is multiplied by 2
    Γ = θ * I(n) + ϵ * Ξ
    Γ11 = @view Γ[1 : l, 1 : l] 
    Γ12 = @view Γ[1 : l, l + 1 : n] 
    Γ22 = @view Γ[l + 1 : n, l + 1 : n]
    α = 1 / ϵ * maximum(eigvals(Γ11 - Γ12 * inv(Γ22) * Γ12')) * 2   # The threshold is multiplied by 2
    A = zeros(n, n)
    A[diagind(A)[1 : l]] .= ϵ * α
    E = ϵ * Ξ
    G = [zeros(n + 1) vcat(zeros(1, n), E)]
    H = [zeros(1, n + 1); diagm(diag(A))' * [-ones(n) diagm(ones(n))]]
    Network(node, G - H, P, [rand(d)*1e-3; x0])
end 

function undirectedcluster(node::AbstractNodeDynamics, 
                           Φ::AbstractMatrix, 
                           γ::AbstractVector, 
                           Δ::AbstractMatrix, 
                           P::AbstractMatrix,
                           J::AbstractVector, 
                           x0::AbstractVector=rand(size(Φ,1)*size(P,1))*1e-3)
    l = length(J) 
    d = size(P, 1) 
    n = size(Φ, 1) 
    Gs = clusters(J)
    ks = length.(Gs)
    (typeof(Φ) <: BlockMatrix) || (Φ = BlockArray(Φ, ks, ks))
    num = maximum(diag(Δ)) + 2 * (l - 1) * maximum([μ(getblock(Φ, i, j)) for i in 1 : l, j in 1 : l if i ≠ j])
    denum = map(1 : l) do i 
            mat = copy(getblock(Φ, i, i))
            mat[end] -= γ[i]
            -maximum(eigvals(mat))
        end |> maximum
    β = num / denum 
    E = Φ
    foreach(i ->  setblock!(E, β * getblock(E, i, i), i, i), 1 : l) 
    α = β * γ
    E = [zeros(l, l) zeros(l, n); zeros(n, l) E]
    H1 = zeros(n, l)
    for (k, (j, αi)) in zip(cumsum(ks), enumerate(α))
        H1[k, j] = -αi
    end
    H2 = zeros(n, n)
    for (k, αi) in zip(cumsum(ks), α)
        H2[k, k] = αi
    end
    H = [zeros(l, l + n); [H1 H2]]
    Network(node, E - H, P, [rand(d*l)*1e-3; x0])
end 

μ(Ψ::AbstractMatrix) = 1 / 2 * maximum(size(Ψ)) * maximum(abs.(Ψ)) 

function clusters(J::AbstractVector)
    J = sort(J)
    i = 1 
    map(J) do j 
        indices = i : j
        i = j + 1 
        indices 
    end
end
