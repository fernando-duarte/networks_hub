using Pkg

Pkg.add(Pkg.PackageSpec(;name="Distributions", version="0.24.15"))
Pkg.add("Quadrature")

using Quadrature, Distributions, Random

# Setting seed
Random.seed!(123)

## Number of nodes
N = 3
# N = prase(Int64, ARGS[1]) - if called via bash 

#User Inputs
α0 = 10 .* (1 .- rand(N))
reltol_val = 1e-10
abstol_val = 1e-10
alg = HCubatureJL()
batch_val = 0
parallel_val = 0
time_start = time()
mem_start = Sys.free_memory()/1000/2^20

# Mean
dist_dirichlet_pdf(x,p) = exp(Distributions.logpdf(Distributions.Dirichlet(p),x)) 
f_dirichlet(x,p) = dist_dirichlet_pdf([x;1.00-sum(x)],p) .* [x;1.00-sum(x)]
prob = QuadratureProblem(f_dirichlet,zeros(N-1),ones(N-1), α0)
#prob = QuadratureProblem(f, zeros(N) , [1000.0, 1000.0], α0)
sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)

mean_dirichlet(p) = p./sum(p)

mean_dirichlet(α0)
