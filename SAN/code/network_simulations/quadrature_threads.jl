import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using Quadrature, Cuba, Cubature, Base.Threads
using Distributions, Random
using DataFrames, CSV
using Flux, CUDA

## User Inputs
N = 6
tol = 5
ind_gpu = 0 # indicator whether to use gpu
alg = CubaSUAVE() # CubaDivonne() #works for CubaCuhre, CubaDivonne CubaSUAVE, fails for cubavega

display("N = $N")
# Setting up Variables
if ind_gpu == 1
    α0 = 10 .* (1 .- rand(N)) |> gpu
else
    α0 = 10 .* (1 .- rand(N))
end
reltol_val = 10.0^(-tol)
abstol_val = 10.0^(-tol)

# Setting up function
dist_dirichlet_pdf(x,p) = Distributions.pdf(Distributions.Dirichlet(p),x)

#=
function f_dirichlet(dx,x,p)
    Threads.@threads for i in 1:N
        dx[i] = (dist_dirichlet_pdf([x;1.00-sum(x)],p) .* [x;1.00-sum(x)])[i]
    end
end

=#
function f_dirichlet(x,p)
    sum((dist_dirichlet_pdf([x;1.00-sum(x)],p) .* [x;1.00-sum(x)]))
end


function f_dirichlet_threads(dx,x,p)
    w = [x;1.00-sum(x)]
  Threads.@threads for i in 1:size(w,2)
   # dx[i] = sum(sin.(@view(x[:,i])))
    dx[i] = sum((dist_dirichlet_pdf(@view(w[:,i]),p) .* @view(w[:,i])))
  end
end

# Solving Integral
prob = QuadratureProblem(f_dirichlet,zeros(N-1),ones(N-1), α0, nout = 1)
prob2 = QuadratureProblem(f_dirichlet_threads,zeros(N-1),ones(N-1), α0, nout = 1) 
time_start = time()
mem_usage = @allocated @time sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)
mem_usage = @allocated @time sol_mean2 = Quadrature.solve(prob2,alg,reltol=reltol_val,abstol=abstol_val)

total_mem = mem_usage/1000/2^20

# Checking Answer
display(sol_mean.u)
display(sol_mean2.u[1])
test_passed = (sol_mean.u - 1.00 < 0.001) + (sol_mean2.u[1] - 1.00 < 0.001)

display(test_passed)

## Saving Meta Data
time_end = time()
total_time = time_end-time_start
meta_data = DataFrame(nodes = N, tolerance = abstol_val, gpu = ind_gpu, total_time_min = total_time/60, total_mem_gb = total_mem, tests_passed = test_passed, total_tests = N, final_mean = sol_mean[1])
CSV.write("meta_data_$(N)_$(tol).csv", meta_data)
