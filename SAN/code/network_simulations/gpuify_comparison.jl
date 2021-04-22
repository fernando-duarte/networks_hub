import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA

@test Threads.nthreads()>1

M= 10 # number of independent beta random variables
atol=1e-6
rtol=1e-3
nvec=10000000
maxevals=100000000

# Naive GPU comparison
function int(x, f)
   f[1] = pdf(Product(Beta.(1.0,2.0*ones(M))),x)
end

function int_gpu(x, f)
   CUDA.@sync f[1] = pdf(Product(Beta.(1.0,(2.0*ones(M)))),CuArray(x)) 
end

# gpu appears slower
display(@benchmark cuhre(int, M, 1, atol=atol, rtol=rtol))
display(@benchmark cuhre(int_gpu, M, 1, atol=atol, rtol=rtol))
display(@benchmark divonne(int, M, 1, atol=atol, rtol=rtol,maxevals=maxevals)[1])
display(@benchmark divonne(int_gpu, M, 1, atol=atol, rtol=rtol,maxevals=maxevals)[1])

# Now Comparing with  FLoops.jl
function int_thread_col(x, f)
     Threads.@threads for i in 1:size(x,2)
      f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),@view(x[:,i]))
    end
end

function int_thread_col_floop(x, f)
    @floop for i in 1:size(x,2)
      f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),@view(x[:,i]))
    end
end

function int_thread_col_floop_cpu(x, f, ex = ThreadedEx())
    @floop ex for i in 1:size(x,2)
        f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),@view(x[:,i]))
    end
end

function int_thread_col_floop_gpu(x, f, ex = has_cuda_gpu() ? CUDAEx() : ThreadedEx())
    @floop for i in 1:size(x,2)
        f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),CuArray(x[:,i]))
    end
end

display(@benchmark cuhre(int_thread_col, M, 1, atol=atol, rtol=rtol))
display(@benchmark cuhre(int_thread_col_floop, M, 1, atol=atol, rtol=rtol)) #floop seems faster than Threads.@threads
display(@benchmark cuhre(int_thread_col_floop_cpu, M, 1, atol=atol, rtol=rtol)) 
display(@benchmark cuhre(int_thread_col_floop_gpu, M, 1, atol=atol, rtol=rtol)) #floop and gpu seem in between in terms of speed 
