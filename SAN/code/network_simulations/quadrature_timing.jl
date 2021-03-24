using Pkg
Pkg.add("Quadrature")
Pkg.add("DistributionsAD")
Pkg.add("SpecialFunctions")
Pkg.add("ForwardDiff")
Pkg.add("FiniteDiff")
Pkg.add("Zygote")
Pkg.add("Cuba")
Pkg.add("Cubature")
Pkg.add("HCubature")
Pkg.add("CSV")
Pkg.add("DataFrames")

using Quadrature, DataFrames, CSV
using LinearAlgebra, Random, Distributions, DistributionsAD
using ForwardDiff, FiniteDiff, Zygote, Cuba, Cubature, HCubature
using Distributed

# Setting seed
Random.seed!(123)

## Number of nodes
N =2 
# N = prase(Int64, ARGS[1]) - if called via bash 

#User Inputs
α0 = 10 .* (1 .- rand(N))
β0 = 50 .* (1 .- rand(N))
reltol_val = 1e-5
abstol_val = 1e-5
alg = HCubatureJL()
batch_val = 0
parallel_val = 0
time_start = time()
mem_start = Sys.free_memory()/1000/2^20

## Dirichlet distribution
dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.TuringDirichlet(N,p[1]),x[1]))
f(x,p) = dist_pdf(x,p) .* x
prob = QuadratureProblem(f,0.0,1.0,[α0[1]...])
sol = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)[1]

# Matrix Beta
dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.MatrixBeta(N,p[1],p[2]),x[1]))
f(x,p) = dist_pdf(x,p) .* x
prob = QuadratureProblem(f,0.0,1.0,[α0[1]..., β0[1]...])
sol = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)[1]

## Mean of the Sum (equivalent to sum of the mean)
dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))
f(x,p) = dist_pdf(x,p) .* x
final_sum = 0

prob = QuadratureProblem(f,0.0,1.0,[α0[1]...,β0[1]...])
sol = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)[1]
    

for i = 1:N
    prob = QuadratureProblem(f,0.0,1.0,[α0[i]...,β0[i]...])
    sol = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)[1]
    final_sum += sol
end
    
## Sum of the Variance
dist_pdf(x,p) = DistributionsAD.pdf(DistributionsAD.Beta(p[1],p[2]),x[1])
final_var = 0
mean_val(α0, β0, j) = α0[j]/(α0[j] + β0[j])


for j = 1:N
    f(x,p) = dist_pdf(x,p) .* (x .- mean_val(α0, β0, j))^2
    prob = QuadratureProblem(f,0.0,1.0,[α0[j]...,β0[j]...])
    sol = Quadrature.solve(prob,alg,reltol=1e-5,abstol=1e-5)[1]
    final_var += sol
end

## Sum of the Squarteroot of the square. 
dist_pdf(x,p) = DistributionsAD.pdf(DistributionsAD.Beta(p[1],p[2]),x[1])
f(x,p) = (dist_pdf(x,p) .* (x^2))^1/2
final_sqrt = 0
for i = 1:N
    prob = QuadratureProblem(f,0.0,100.0,[α0[i]...,β0[i]...])
    sol = Quadrature.solve(prob,alg,reltol=1e-5,abstol=1e-5)[1] # Expected value times number of observations is the mean of the sum
    final_sqrt += sol
end

time_end = time()
mem_end = Sys.free_memory()/1000/2^20
total_time = time_end - time_start
total_mem = mem_start - mem_end

data_out = DataFrame()
data_out.N = [N]
data_out.algorithm = string(alg)
data_out.reltol = reltol_val
data_out.abstol = abstol_val
data_out.batch = batch_val
data_out.parallelism = parallel_val
data_out.computation_time_sec = total_time
data_out.maxRAM_used_gb = total_mem

CSV.write("quadrature_results.csv",data_out)