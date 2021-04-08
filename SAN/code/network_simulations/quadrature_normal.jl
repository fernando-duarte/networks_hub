import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using Quadrature, Cuba, Cubature, Base.Threads
using Distributions, Random, LinearAlgebra
using DataFrames, CSV

#User Input
N = 8
tol = 3
alg = HCubatureJL()

# Setting up problem
reltol_val = 10.0^(-tol)
abstol_val = 10.0^(-tol)

means = rand(N)
vars = Diagonal(abs.(rand(N,N)))
d = MvNormal(means , vars)
m(x , p) = pdf(d, x) .* x
prob = QuadratureProblem(m, -Inf*ones(N) , Inf*ones(N))
time_start = time()
mem_usage = @allocated sol_mean = Quadrature.solve(prob,alg,reltol=1e-2,abstol=1e-2)
total_mem = mem_usage/1000/2^20

display(sol_mean)
display(means)


test_passed = sum((abs.(sol_mean .- means) .< 1e-2))

## Saving Meta Data
time_end = time()
total_time = time_end-time_start
meta_data = DataFrame(nodes = N, tolerance = abstol_val, total_time_min = total_time/60, total_mem_gb = total_mem, tests_passed = test_passed, total_tests = N, final_mean = sol_mean[1], final_mean_correct = means[1])
CSV.write("meta_data_$(N)_$(tol).csv", meta_data)
