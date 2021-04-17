using Convex, SCS, LinearAlgebra

# Construct model 
A = [0 1; -2 -3]
n = size(A, 1)
Q = I(n)
n = size(A, 1)
P = Semidefinite(n)
probconstriants = [
    A' * P + P * A ≤ -Q, 
    A' * P + P * A ≥ -Q
] 
problem = satisfy(probconstriants)

# Optimize model 
solve!(problem, () -> SCS.Optimizer(verbose=false))

# Test 
@show P.value
@show A' * P.value + P.value * A ≈ -Q;
