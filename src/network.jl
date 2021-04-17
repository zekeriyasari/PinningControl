
struct Network{T1<:AbstractNodeDynamics, T2<:AbstractMatrix, T3<:AbstractMatrix}
    node::T1
    outermat::T2
    innermat::T3
    state::Vector{Float64}
end 
Network(node::AbstractNodeDynamics, outermat::AbstractMatrix, innermat::AbstractMatrix) = 
    Network(node, outermat, innermat, rand(size(outermat, 1) * size(innermat, 1)))

function netrhs(dx, x, (f, E, P), t)
    for idx in Iterators.partition(1 : length(x), size(P, 1))
        f(view(dx, idx), view(x, idx), nothing, t)
    end 
    dx .+= kron(E, P) * x
end 

netprob(net::Network, ti::Real, tf::Real) = ODEProblem(netrhs, net.state, (ti, tf), (net.node, net.outermat, net.innermat))
netsolve(net::Network, ti::Real, dt::Real, tf::Real) = solve(netprob(net, ti, tf), saveat=dt)
