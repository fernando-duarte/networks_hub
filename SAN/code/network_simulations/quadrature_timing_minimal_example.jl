# RUN FIRST TIME ONLY
#=
using Pkg
Pkg.generate("quadrature_timing") # - run first time
Pkg.activate("quadrature_timing")
Pkg.add("Quadrature")
Pkg.add(Pkg.PackageSpec(;name="Distributions", version="0.24.15"))
Pkg.add("Random")
Pkg.add("DataFrames")
Pkg.add("CSV")
=#

using Pkg
Pkg.activate("quadrature_timing")
using Quadrature, Distributions, Random, DataFrames, CSV

# Setting seed
Random.seed!(123)

## Number of nodes
N = 6
# N = parse(Int64, ARGS[1]) - if called via bash 

#User Inputs
α0 = 10 .* (1 .- rand(N))
reltol_val = 1e-5
abstol_val = 1e-5
alg = HCubatureJL()
batch_val = 0
parallel_val = 0
time_start = time()
mem_start = Sys.free_memory()/1000/2^20

# Mean
dist_dirichlet_pdf(x,p) = Distributions.pdf(Distributions.Dirichlet(p),x)
f_dirichlet(x,p) = dist_dirichlet_pdf([x;1.00-sum(x)],p) .* [x;1.00-sum(x)]
prob = QuadratureProblem(f_dirichlet,zeros(N-1),ones(N-1), α0)
sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

mean_dirichlet(p) = p./sum(p)

mean_dirichlet(α0)

# Variance
dist_dirichlet_pdf(x,p) = Distributions.pdf(Distributions.Dirichlet(p),x)
f(x,p) = dist_dirichlet_pdf([x;1.00-sum(x)],p) * (([x;1.00-sum(x)] .- mean_dirichlet(α0))).^2
prob = QuadratureProblem(f,zeros(N-1),ones(N-1), α0)
sol_var = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

var_dirichlet(p) = ((p./sum(p)) .- (p./sum(p)).^2) ./(sum(p)+1)
var_dirichlet(α0)

covar_dirichlet(p, i, j) = -p[i]*p[j]/(sum(p)^3 + sum(p))   #only when i != j
covar_dirichlet(α0,2,1)

test_passed = sum((sol_mean .- mean_dirichlet(α0) .< 1e-5) + (sol_var .- var_dirichlet(α0) .< 1e-5))

## Saving Meta Data
mem_end = Sys.free_memory()/1000/2^20
total_mem = mem_start - mem_end 
time_end = time()
total_time = time_end-time_start
meta_data = DataFrame(nodes = N, total_time_min = total_time/60, total_mem_gb = total_mem, tests_passed = test_passed, total_tests = N*2)
CSV.write("meta_data_$N.csv", meta_data)
