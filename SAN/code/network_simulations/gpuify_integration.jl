import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA
using SpecialFunctions

@test Threads.nthreads()>1

# User Inputs
M= 25 # number of independent uniform random variables
atol=1e-6
rtol=1e-3
nvec=1000000
maxevals=100000000

# Initializing Matrices
ones_mat = CuArray(ones(Float32, M))
result = CUDA.ones(Float32, (nvec,1))
x_cpu = rand(Float64, (nvec, M))
x = CuArray(rand(Float32, (nvec, M)))
x_1d_gpu = CuArray(rand(Float64, M))
x_1d_cpu = rand(Float64, M)

# Initializing Functions
function int_cpu(x, f)
   f[1] = pdf(Product(Beta.(1.0,2.0*ones(M))),x)
end

function int_cpu2(x, f)
   f[1] = vec(prod(x'.^(1.0-1.0) .* (1.0 .- x').^(2.0-1.0)./(gamma(1.0)*gamma(2.0)/gamma(3.0)),dims=2))[1]
end

function int_thread_col_cpu(x, f)
    Threads.@threads for i in 1:size(x,2)
      f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),@view(x[:,i]))
    end
end

function int_thread_col_cpu2(x, f)
    @floop for i in 1:size(x,2)
      f[i] = prod(@view(x[:,i]).^(1.0f0-1.0f0) .* (1.0f0 .- @view(x[:,i])).^(2.0f0-1.0f0)./(gamma(1.0f0)*gamma(2.0f0)/gamma(3.0f0)),dims=1)[1]
    end
end

function int_thread_el_cpu(x,f)
   f[1,:] .= 1.0
   Threads.@threads for j in 1:size(x,2)
       for i in 1:size(x, 1)
           f[1, j] *= pdf(Beta(1.0,2.0),@view(x[i,j]))
       end
   end
end

function int_thread_el_cpu2(x,f)
    f[1,:] .= 1.0 
    power1 = 1.0f0 - 1.0f0
    power2 = 2.0f0 - 1.0f0
    denominator = (exp.(CUDA.lgamma.(1.0f0)) .* exp.(CUDA.lgamma.(2.0f0)) ./ exp.(CUDA.lgamma.(3.0f0)))[1]
    Threads.@threads for j in 1:size(x,2)
       @floop for i in 1:size(x, 1)
           f[1, j] *= ((x[i,j])^power1 * (1.0f0 - x[i,j])^power2) /denominator
       end
   end
end

# Benchmarking CPU options

# cuhre
display(@benchmark cuhre($int_cpu, $M, 1, atol=$atol, rtol=$rtol)) # (2.0 ms for M = 5, 650ms for M=15)
display(@benchmark cuhre($int_cpu2, $M, 1, atol=$atol, rtol=$rtol)) # (500 mus for M = 5, 100ms for M = 15, 38s for M = 25)

#suave 
display(@benchmark suave($int_thread_col_cpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 60 ms for M = 5, nvec = 1000000, 35s M = 15
display(@benchmark suave($int_thread_col_cpu2, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #15 ms for M = 5, nvec = 1000000, 7s M = 15, 27ms for M = 25 but gets incorrect answer
display(@benchmark suave($int_thread_el_cpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 45 ms for M = 5, nvec = 1000000, 21s M = 15
display(@benchmark suave($int_thread_el_cpu2, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 20ms for M = 5, nvec = 1000000, 8s M = 15

#divonne
display(@benchmark divonne($int_thread_col_cpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 380 ms for M = 5, nvec = 1000000
display(@benchmark divonne($int_thread_col_cpu2, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 4.5 ms for M = 5, nvec = 1000000
display(@benchmark divonne($int_thread_el_cpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 440 ms for M = 5, nvec = 1000000
display(@benchmark divonne($int_thread_el_cpu2, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 430 ms for M = 5, nvec = 1000000

=#
## Takeaway: int_col and int_el work well for suave and divonne, but int_cpu works best for cuhre. 
## Takeaway2: Writing out the formula is much faster than using the product. 

## Integration with GPU
function beta_pdf_gpu(x, a, b)
     prod(x.^(a-1.0f0) .* (1.0f0 .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function int_gpu(x, f)
   f[1] = vec(beta_pdf_gpu(CuArray(x),1.0f0,2.0f0))[1]
end

function int_thread_col_gpu(x, f)
    @floop for i in 1:size(x,2)
      f[i] = vec(beta_pdf_gpu(CuArray(x[:,i]),1.0f0,2.0f0))[1]
    end
end

function int_thread_col_gpu_precompile(x, f)
    power1 = 1.0f0 - 1.0f0
    power2 = 2.0f0 - 1.0f0
    denominator = (exp.(CUDA.lgamma.(1.0f0)) .* exp.(CUDA.lgamma.(2.0f0)) ./ exp.(CUDA.lgamma.(3.0f0)))[1]
    ones_mat = ones(Float32, (size(x)))
    val2 = (x.^(power1) .* (ones_mat .- x).^(power2))./denominator
    @floop for i in 1:size(x,2)
      f[i] = reduce(*,@view(val2[:,i]))
    end
    return f
end

## Benchmarking GPU Options

# cuhre
display(@benchmark cuhre($int_gpu, $M, 1, atol=$atol, rtol=$rtol)) # 70 ms for M = 5, 11.7 s for M = 15)

#suave
display(@benchmark suave($int_thread_col_gpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #2.3 s for M = 5, nvec = 1000000, way too long for M = 15
display(@benchmark suave($int_thread_col_gpu_precompile, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #12 ms  for M = 5, nvec = 10000000. 5.8s for M = 15


#divonne
display(@benchmark divonne($int_thread_col_gpu, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #1.5s for M = 5, nvec = 1000000, way too long for M = 15
display(@benchmark divonne($int_thread_col_gpu_precompile, $M, 1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #3 ms  for M = 5, nvec = 10000000, 17s for M = 15

