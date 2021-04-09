import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Distributed
addprocs(2)
@show nprocs(), workers()

@everywhere module TestQuadrature
    
    using Quadrature, Cuba, Cubature, Base.Threads
    using Distributions, Random
    using DataFrames, CSV
    using Distributed
    
    N = 3
    tol =5
    ind_gpu = 0 # indicator whether to use gpu
    alg = CubaSUAVE() # CubaDivonne() #works for CubaCuhre, CubaDivonne CubaSUAVE, fails for cubavega

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
        dx[i] = sum((dist_dirichlet_pdf(@view(w[:,i]),p) .* @view(w[:,i])))
      end
    end

    function f_dirichlet_workers(dx,x,p)
        w = [x;1.00-sum(x)]
        @sync @distributed for i in 1:size(w,2)
        dx[i] = sum((dist_dirichlet_pdf(@view(w[:,i]),p) .* @view(w[:,i])))
      end
    end

    function run()
    # Solving Integral
        prob = QuadratureProblem(f_dirichlet,zeros(N-1),ones(N-1), α0, nout = 1)
        prob2 = QuadratureProblem(f_dirichlet_threads,zeros(N-1),ones(N-1), α0, nout = 1) 
        prob3 = QuadratureProblem(f_dirichlet_workers,zeros(N-1),ones(N-1), α0, nout = 1) 
        time_start = time()
        mem_usage = @allocated @time sol_mean = Quadrature.solve(prob,alg,reltol=reltol_val,abstol=abstol_val)
        mem_usage = @allocated @time sol_mean2 = Quadrature.solve(prob2,alg,reltol=reltol_val,abstol=abstol_val)
        mem_usage = @allocated @time sol_mean3 = Quadrature.solve(prob3,alg,reltol=reltol_val,abstol=abstol_val)
    display(sol_mean)
    display(sol_mean2)
    display(sol_mean3)
    end
end

using .TestQuadrature
TestQuadrature.run()
