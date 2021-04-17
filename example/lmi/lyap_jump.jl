using JuMP, SCS

# Construct model 
A = [0 1; -2 -3] 
n = size(A, 1)
Q = collect(I(n))
model = Model(SCS.Optimizer)
set_silent(model)

# Add contraints 
@variable(model, P[1:n, 1:n], PSD)
@SDconstraint(model, P ≥ 0)
@SDconstraint(model, A'*P + P*A ≥ -Q)
@SDconstraint(model, A'*P + P*A ≤ -Q)

# Optimize model 
optimize!(model)

# Test 
P = value.(P)
@show P
@show A'*P + P*A ≈ -Q;
