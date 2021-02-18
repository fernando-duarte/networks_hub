using Pkg
Pkg.add(Pkg.PackageSpec(;name="InfiniteOpt", version="0.3.2"))

using LinearAlgebra, DataFrames, XLSX, JuMP, Ipopt, Random, Distributions
using Plots, InfiniteOpt
plotly(ticks=:native)  

N = 3 # keep largest `N' nodes by assets

c = [1.95, 1.75, 1.65]
b = [0.6, 9.7, 0.9]
p_bar = [1.0, 0.9, 0.9]
assets = [2.0, 1.8, 1.7]
delta = [0.1, 0.05, 0.02]
w= [0.16, 0.15, 0.14]
g0 = 0.05   # bankruptcy cost

# initial guess
rng =  Random.seed!(123)
A0 = zeros(N,N)
c0 =  c
b0 = b
α0 = 1.0*ones(N)
β0= 50.0*ones(N) 

## Setting up model
m = InfiniteModel(Ipopt.Optimizer);

# Parameter everything depends upon (shocks to each node)
@infinite_parameter(m, x[i=1:N] in Beta(α0[i],β0[i]), num_supports = 1000) #parameter of shocks

# setting up variables, starting position, and constraints
@infinite_variable(m, p[i=1:N](x), start = 0) # clearing vector. 
@hold_variable(m, c[i=1:N] >= 0) #outside assets for node i, 
@hold_variable(m, b[i=1:N] >= 0) #outside liabilities for node i
@hold_variable(m, A[i=1:N, j=1:N]>= 0) #net payments scaled by total liabilites for firm i between nodes i and j
@hold_variable(m, α[i=1:N]>= 0) # parameter that determines the shocks
@hold_variable(m, β[i=1:N]>= 0) # parameter that determines shocks 

for i = 1:N
    @constraint(m, w[i] <= c[i] <= assets[i]) # outside assets greater than networth and lower than total assets
    @constraint(m, 0 <= b[i] <= p_bar[i]) #outside liabilities greater than 0 and less than total liabilities
    @constraint(m, 0.8 <= α[i] <= 3.0)
    @constraint(m, 1.0 <= β[i] <= 150.0)
end

@constraint(m, -sum(A,dims=2).*p_bar .+ p_bar .==  b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * p_bar .== assets .- c ) # payments from other nodes add up to inside assets d

for i = 1:N
    for j = 1:N
        @constraint(m, 0<= A[i, j] <= 1)
    end
end

for i = 1:N
    @constraint(m, A[i,i]==0)
end

for i=1:N
    j=1
    while j < i
        @constraint(m, 0 == A[i,j] * A[j,i])
        j += 1
    end
end

# Clearing Vector Constraints
@constraint(m, p .== transpose(A)*p + c - x)

# Setting up objective function
@objective(m, Min, expect(sum(c.*x + p_bar - p),x)) # min/max E[(ci * xi + p_bar i - pi(x)))]

# Training Model
@time JuMP.optimize!(m);

## Checking Solution 
# variable values after training
csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)
tol = 1e-5

#tests 
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(-tol .<=Asol0.<=1+tol)
    @test all(-tol .<=bsol0.<=data.p_bar)
    @test all(-tol .<=csol0.<=data.assets)   
end
