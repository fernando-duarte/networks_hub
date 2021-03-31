using Pkg
Pkg.activate("quadrature_timing")
Pkg.instantiate
using Quadrature, DataFrames, CSV
using LinearAlgebra, Random, Distributions, DistributionsAD
using ForwardDiff, FiniteDiff, Zygote, Cuba, Cubature, HCubature
using Distributed
using Test
using DistributionsAD: TuringDirichlet

# Setting seed
Random.seed!(123)

## Number of nodes
N = 2
# N = prase(Int64, ARGS[1]) - if called via bash 

#User Inputs
α0 = 10 .* (1 .- rand(N))
β0 = 50 .* (1 .- rand(N))
reltol_val = 1e-3
abstol_val = 1e-3
alg = HCubatureJL()
batch_val = 0
parallel_val = 0
time_start = time()
mem_start = Sys.free_memory()/1000/2^20

## Dirichlet distribution
function dist_dirichlet_pdf_corrected(x,p)
    if sum(x) - 1.00 < 1e6
        exp(DistributionsAD.logpdf(DistributionsAD.TuringDirichlet(p),x)) 
    else
        0.0
    end
end

# Mean
dist_dirichlet_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.TuringDirichlet(p),x)) 
f(x,p) = dist_dirichlet_pdf(x,p) * sum(x)
prob = QuadratureProblem(f,zeros(N),ones(N), α0)
#prob = QuadratureProblem(f, zeros(N) , [1000.0, 1000.0], α0)
sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

mean_dirichlet(p) = p./sum(p)

mean_dirichlet(α0)

# Variance
dist_dirichlet_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.TuringDirichlet(p),x)) 
f(x,p) = dist_dirichlet_pdf(x,p) * sum(((x .- mean_dirichlet(α0))).^2)
prob = QuadratureProblem(f,zeros(N),ones(N), α0)
#prob = QuadratureProblem(f, zeros(N) , [1000.0, 1000.0], α0)
sol_var = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

var_dirichlet(p) = ((p./sum(p)) .- (p./sum(p)).^2) ./(sum(p)+1)
var_dirichlet(α0)

covar_dirichlet(p, i, j) = -p[i]*p[j]/(sum(p)^3 + sum(p))   #only when i != j
covar_dirichlet(α0,2,1)

# Mean of Sqrt
dist_dirichlet_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.TuringDirichlet(p),x)) 
f(x,p) = dist_dirichlet_pdf(x,p) * sum(sqrt.(x))^2
prob = QuadratureProblem(f,zeros(N),ones(N), α0)
#prob = QuadratureProblem(f, zeros(N) , [1000.0, 1000.0], α0)
sol_sqrt = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

#=
# Matrix Beta

dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.MatrixBeta(N,p[1],p[2]),x))
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

## Tests
# is this broken? 
w = rand(DistributionsAD.TuringDirichlet(α0),5) 

@testset "TuringDirichlet" begin
        dim = 3
        n = 4
        for alpha in (2, rand())
            d1 = TuringDirichlet(dim, alpha)
            d2 = Dirichlet(dim, alpha)
            #d3 = TuringDirichlet(d2)
            @test d1.alpha == d2.alpha # == d3.alpha
            @test d1.alpha0 == d2.alpha0 # == d3.alpha0
            @test d1.lmnB == d2.lmnB # == d3.lmnB

            s1 = rand(d1)
            @test s1 isa Vector{Float64}
            @test length(s1) == dim

            s2 = rand(d1, n)
            @test s2 isa Matrix{Float64}
            @test size(s2) == (dim, n)
        end

        for alpha in (ones(Int, dim), rand(dim))
            d1 = TuringDirichlet(alpha)
            d2 = Dirichlet(alpha)
            #d3 = TuringDirichlet(d2)
            @test d1.alpha == d2.alpha # == d3.alpha
            @test d1.alpha0 == d2.alpha0 # == d3.alpha0
            @test d1.lmnB == d2.lmnB # == d3.lmnB

            s1 = rand(d1)
            @test s1 isa Vector{Float64}
            @test length(s1) == dim

            s2 = rand(d1, n)
            @test s2 isa Matrix{Float64}
            @test size(s2) == (dim, n)
        end
end


## Turing Dirichlet 
DistributionsAD.logpdf(TuringDirichlet([1.2,3.4,5.1]),[0.1,0.2,0.5])

Zygote.gradient( (alpha1,alpha2,alpha3) -> exp(DistributionsAD.logpdf(TuringDirichlet([alpha1,alpha2,alpha3]),[0.1,0.2,0.5])) , (1.2,1.3,1.5)...)

zg(a1,a2,x1,x2) = Zygote.gradient( (alpha1,alpha2)-> exp(DistributionsAD.logpdf(TuringDirichlet([alpha1,alpha2]),[x1,x2])) , (a1,a2)...)
zg(1.3,1.5,0.45,0.12)

using SpecialFunctions
pdfDirichlet(a1,a2,x1,x2) = exp(DistributionsAD.logpdf(TuringDirichlet([a1,a2]),[x1,x2]))
mg(a1,a2,x1,x2)=pdfDirichlet(a1,a2,x1,x2)*(digamma(a1+a2)-digamma(a1)+log(x1))
zg(1.3,1.5,0.45,0.12)
mg(1.3,1.5,0.45,0.12)

=#
