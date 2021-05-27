import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using CUDA, Test, BenchmarkTools
using SpecialFunctions

## Initializing Matricies
M = 15
nvec = 5000
x_cpu = rand(Float64, (M, nvec))
x = convert(CuArray{Float32},x_cpu)
result = CUDA.ones(Float32, (M,nvec))

function beta_pdf_gpu_basic!(x,a,b,denom,result)
    local tid = 1
    while (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))./denom
        tid += 1
    end
    return nothing
end

function beta_pdf_gpu_basic(x, a, b, result)
    denom = gamma(a)*gamma(b)/gamma(a+b)
    @cuda beta_pdf_gpu_basic!(x, a, b,denom ,result)
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu_blocks!(x,a,b,denom,result)
    tid = blockIdx().x
    if (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))./denom
    end
    return nothing
end

function beta_pdf_gpu_blocks(x, a, b, result,numblocks)
    denom = gamma(a)*gamma(b)/gamma(a+b)
    @cuda blocks=numblocks beta_pdf_gpu_blocks!(x, a, b,denom ,result)
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu_threads!(x,a,b,denom,result)
    tid = threadIdx().x
    if (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))./denom
    end
    return nothing
end

function beta_pdf_gpu_threads(x, a, b, result,numthreads)
    denom = gamma(a)*gamma(b)/gamma(a+b)
    @cuda threads=numthreads beta_pdf_gpu_threads!(x, a, b,denom ,result)
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu(x::CuArray,a,b,dim)
    prod(x.^(a-1.0f0) .* (1.0f0 .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=dim)
end

@btime beta_pdf_gpu($x, $1.0f0, $2.0f0,$1)
@btime beta_pdf_gpu_basic($x, $1.0f0, $2.0f0, $result)
@btime beta_pdf_gpu_blocks($x, $1.0f0, $2.0f0, $result, $10)
@btime beta_pdf_gpu_threads($x, $1.0f0, $2.0f0, $result, $100)

# For M = 15 and nvec = 50 the new versions are faster than the old way by about twice the speed
# For M = 15 and nvec = 500 new versions are still much faster
# For M = 15 and nvec = 5000 new versions are still much faster 
# For M = 15 and nvec = 50000 new versions are much slower (1s vs 32 mus)
# threads and block version get incorrect answer if nvec too large compared to the number of threads/blocks because the don't do calculation for all columns of x surprisingly






## Attempt to get thread block combination to work but not working yet. 
function beta_pdf_gpu_threads_blocks!(x,a,b,denom,result)
    tid = (threadIdx().x - 1) + (blockIdx().x - 1) * blockDim().x
    if (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))./denom
    end
    return nothing
end

function beta_pdf_gpu_threads_blocks(x, a, b, result,numthreads, numblocks)
    denom = gamma(a)*gamma(b)/gamma(a+b)
    cache = @cuDynamicSharedMem(Int64, threadsPerBlock)
    @cuda threads=numthreads blocks=numblocks beta_pdf_gpu_threads_blocks!(x, a, b,denom ,result)
    final_ans = prod(result,dims=1)
end

@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $result,$5, $2)
