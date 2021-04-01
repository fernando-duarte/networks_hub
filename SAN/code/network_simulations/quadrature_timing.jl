import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Quadrature, Distributions, Random, DataFrames, CSV, Cuba, Cubature

# Setting seed
Random.seed!(123)

## Number of nodes
N = 3
#N = parse(Int64, ARGS[1]) #- if called via bash 

#User Inputs
α0 = 10 .* (1 .- rand(N))
reltol_val = 1e-5
abstol_val = 1e-5
alg = CubaDivonne() #works for CubaCuhre, CubaDivonne CubaSUAVE, fails for cubavegas
batch_val = 0
parallel_val = 0
time_start = time()

mean_dirichlet(p) = p./sum(p)

# Mean
final_mean = zeros(N)

function fun_compute_mean(N)
    for i = 1:N
        dist_dirichlet_pdf(x,p) = Distributions.pdf(Distributions.Dirichlet(p),x)
        f_dirichlet(x,p) = dist_dirichlet_pdf([x;1.00-sum(x)],p) .* [x;1.00-sum(x)][i]
        prob = QuadratureProblem(f_dirichlet,zeros(N-1),ones(N-1), α0,maxiters = 1e9)
        sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)
        global final_mean[i] = sol_mean[1]
    end 
    return final_mean
end 

mem_usage = @allocated fun_compute_mean(N)
total_mem_mean = mem_usage/1000/2^20

# Variance
final_var = zeros(N)
function fun_compute_var(N)
    for i = 1:N
        dist_dirichlet_pdf(x,p) = Distributions.pdf(Distributions.Dirichlet(p),x)
        f(x,p) = dist_dirichlet_pdf([x;1.00-sum(x)],p) * ((([x;1.00-sum(x)] .- mean_dirichlet(α0))).^2)[i]
        prob = QuadratureProblem(f,zeros(N-1),ones(N-1), α0)
        sol_var = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)
        global final_var[i] = sol_var[1]
    end
    return final_var
end 

mem_usage = @allocated fun_compute_var(N)
total_mem_var = mem_usage/1000/2^20
total_mem = total_mem_mean + total_mem_var 

var_dirichlet(p) = ((p./sum(p)) .- (p./sum(p)).^2) ./(sum(p)+1)
var_dirichlet(α0)

covar_dirichlet(p, i, j) = -p[i]*p[j]/(sum(p)^3 + sum(p))   #only when i != j
covar_dirichlet(α0,2,1)

test_passed = sum((abs.(final_mean .- mean_dirichlet(α0)) .< 1e-4) + (abs.(final_var .- var_dirichlet(α0)) .< 1e-4))

## Saving Meta Data
time_end = time()
total_time = time_end-time_start
meta_data = DataFrame(nodes = N, total_time_min = total_time/60, total_mem_gb = total_mem, tests_passed = test_passed, total_tests = N*2)
CSV.write("meta_data_$N.csv", meta_data)

# Solving using cuba directly - works for N = 3 but not other values
#divonne((z,f,ndim=N,nvec=N) -> f[] = Distributions.pdf(Distributions.Dirichlet(α0),[z;1.00-sum(z)]).*[z;1.00-sum(z)][1])


