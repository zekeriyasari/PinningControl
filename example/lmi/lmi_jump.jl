using JuMP, SCS 
using LightGraphs, GraphPlot
using LinearAlgebra

# Construct a graph 
θ = 10
ϵ = 50.
graph = random_regular_graph(10, 4)
# gplot(graph)
Ξ = -collect(laplacian_matrix(graph))

# Construct an optimization model 
model = Model(SCS.Optimizer)
set_silent(model)

# Add variables 
n = size(Ξ, 1)
@variable(model, A[1 : n, 1 : n], PSD)

# Specify the number for pinned nodes 
l = n ÷ 2

# Add constraints 
Γ = θ * I(n) + ϵ * Ξ
σ = 1.
@SDconstraint(model, A + σ * I(n) ≥ 0 )
@SDconstraint(model, Γ - A ≤ 0)
for i in 1 : n 
    for j in 1 : n 
        if i ≠ j 
            @constraint(model, A[i, j] == 0) # Off diagonal elements are zero
        end 
    end 
end 
for i in l + 1 : n
    @constraint(model, A[i, i] == 0) # Diagonal element constraint 
end 

# Optimize the model 
optimize!(model)

# Get the value 
Aval = value.(A)

# Checks 
@show termination_status(model)
@show eigvals(Aval)
@show all(eigvals(Aval) .> 0)
@show all(eigvals(Γ - Aval) .< 0);
scatter(eigvals(Aval))
Aval[abs.(Aval) .≤  7e-6] .= 0.
display(Aval)
