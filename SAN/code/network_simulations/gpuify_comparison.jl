import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA
using SpecialFunctions
using Referenceables: referenceable

@test Threads.nthreads()>1

# User Inputs
M= 5 # number of independent uniform random variables
atol=1e-6
rtol=1e-3
nvec=10000
maxevals=100000000

# Initializing Matrices
x = CuArray(rand(Float32, (M, nvec)))
ones_mat = CuArray(ones(Float32, M))

function beta_pdf_gpu(x, a, b)
     prod(x.^(a-1.0f0) .* (ones_mat .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_reduce(x, a, b)
    @floop for i = 1:size(x,2)
        reduce(*,x[:,i].^(a-1.0f0) .* (ones_mat .- x[:,i]).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
    end
end

function beta_pdf_gpu_array(x, a, b)
   prod(x.^(CuArray([a])-CuArray([1.0f0])) .* (ones_mat .- x).^(CuArray([b])-CuArray([1.0f0]))./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_precompile(x, a, b)
    power1 = a-CuArray([1.0f0])
    power2 = b-CuArray([1.0f0])
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    prod((x.^(power1) .* (ones_mat .- x).^(power2))./denominator,dims=1)
end


function beta_pdf_gpu_precompile_flip(x, a, b)
    power1 = a-CuArray([1.0f0])
    power2 = b-CuArray([1.0f0])
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    ones_mat_flip = ones_mat'
    prod((x.^(power1) .* (ones_mat_flip .- x).^(power2))./denominator,dims=2)
end

function beta_pdf_gpu_precompile_floop(x, a, b) # super slow
    power1 = a-CuArray([1.0f0])
    power2 = b-CuArray([1.0f0])
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    for j in 1:size(x,2)
        @floop for w in referenceable(x[:,j])
            @reduce(xprod *= (w[]^(power1[1]) * (1.0f0 .- w[])^(power2[1]))/denominator)
        end
        f[j] = xprod
    end
    return f
end


function beta_pdf_gpu_precompile_floop2(x, a, b) # super slow
    power1 = a-CuArray([1.0f0])
    power2 = b-CuArray([1.0f0])
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    @floop CUDAEx() for j in 1:size(x,2)
        xprod = prod((x[:,j].^(power1) .* (ones_mat .- x[:,j]).^(power2))./denominator,dims=1)
        f[j] = xprod[1]
    end
    return f
end

# Benchmarking Different Options
display(@benchmark beta_pdf_gpu($x, 1.0f0, 2.0f0))
display(@benchmark beta_pdf_gpu_array($x,1.0f0,2.0f0)) # running on gpu using user written beta function with float 32 and preallocated ones as float32 (faster)
display(@benchmark beta_pdf_gpu_reduce($x,1.0f0,2.0f0)) # slower
display(@benchmark beta_pdf_gpu_precompile($x,CuArray([1.0f0]),CuArray([2.0f0]))) #faster 
display(@benchmark beta_pdf_gpu_precompile_flip($x',CuArray([1.0f0]),CuArray([2.0f0]))) #marginally slower 

f = ones(Float32, (1, size(x,2)))
display(@benchmark beta_pdf_gpu_precompile_floop($x,CuArray([1.0f0]),CuArray([2.0f0]))) 
display(@benchmark beta_pdf_gpu_precompile_floop2($x,CuArray([1.0f0]),CuArray([2.0f0]))) #much slower
display(@benchmark beta_pdf_gpu_ones_32_array_lgamma($cux32,CuArray([1.0f0]),CuArray([2.0f0]))) # running on gpu using user written beta function with float 32 and preallocated ones as float32 and inbounds (same as above) 

## Trying out multiplication strategies  
display(@benchmark reduce(*, x, dims = 1))
display(@benchmark prod(x, dims=1)) #identical to above 

function parallel_multi(x)
    for i in 1:size(x, 2)
        @floop CUDAEx() for w in referenceable(x[:,i])
         @reduce(val *= w[])
        end
        result[i] = val
    end
    return result 
end
result = CUDA.ones(size(x,2),1)
display(@benchmark parallel_multi(x)) #seems way slower
