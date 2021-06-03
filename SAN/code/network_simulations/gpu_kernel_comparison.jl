import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using CUDA, Test, BenchmarkTools
using SpecialFunctions

## Initializing Matricies
const M = 15
const nvec = 10000
const x_cpu = rand(Float64, (M, nvec))
const x = convert(CuArray{Float32},x_cpu)

function beta_pdf_gpu_basic!(x,a,b,denom,result)
    local tid = 1
    while (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))/denom[1]
        tid += 1
    end
    return nothing
end

function beta_pdf_gpu_basic(x, a, b, result,numthreads,numblocks)
     denom = CuArray([gamma(a)*gamma(b)/gamma(a+b)])
    CUDA.@sync begin
        @cuda threads=numthreads blocks=numblocks beta_pdf_gpu_basic!(x, a, b,denom ,result)
    end
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu_blocks!(x,a,b,denom,result)
    tid = blockIdx().x
    if (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))/denom[1]
    end
    return nothing
end

function beta_pdf_gpu_blocks(x, a, b, result,numblocks)
     denom = CuArray([gamma(a)*gamma(b)/gamma(a+b)])
    CUDA.@sync begin
        @cuda blocks=numblocks beta_pdf_gpu_blocks!(x, a, b, denom, result)
    end
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu_threads!(x,a,b,denom,result)
    tid = threadIdx().x
    if (tid <= length(x))
        result[tid] = (x[tid]^(a-1.0f0) * (1.0f0-x[tid])^(b-1.0f0))/denom[1]
    end
    return nothing
end

function beta_pdf_gpu_threads(x, a, b, result,numthreads)
     denom = CuArray([gamma(a)*gamma(b)/gamma(a+b)])
    @cuda threads=numthreads beta_pdf_gpu_threads!(x, a, b,denom ,result)
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu_threads_blocks!(x,a,b,denom,result)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(x)
        @inbounds result[i] = (x[i]^(a-1.0f0) * (1.0f0-x[i])^(b-1.0f0))/denom[1]
    end
    return nothing
end

function beta_pdf_gpu_threads_blocks(x, a, b, result, threadsPerBlock, blocksPerGrid)
    denom = CuArray([gamma(a)*gamma(b)/gamma(a+b)])
    CUDA.@sync begin
        @cuda threads = threadsPerBlock blocks = blocksPerGrid beta_pdf_gpu_threads_blocks!(x, a, b, denom ,result)
    end
    final_ans = prod(result,dims=1)
end

function beta_pdf_gpu(x::CuArray,a,b,dim)
    prod(x.^(a-1.0f0) .* (1.0f0 .- x).^(b-1.0f0)/(gamma(a)*gamma(b)/gamma(a+b)),dims=dim)
end

# Testing methods are equivalent 
@test beta_pdf_gpu(x, 1.0f0, 2.0f0, 1) ≈ beta_pdf_gpu_basic(x, 1.0f0, 2.0f0, CUDA.ones(Float32, (M,nvec)), 1, 1)
@test beta_pdf_gpu(x, 1.0f0, 2.0f0, 1) ≈ beta_pdf_gpu_blocks(x, 1.0f0, 2.0f0, CUDA.ones(Float32, (M,nvec)), 10)
@test beta_pdf_gpu(x, 1.0f0, 2.0f0, 1) ≈ beta_pdf_gpu_threads(x, 1.0f0, 2.0f0, CUDA.ones(Float32, (M,nvec)), length(x))
@test beta_pdf_gpu(x, 1.0f0, 2.0f0, 1) ≈ beta_pdf_gpu_threads_blocks(x, 1.0f0, 2.0f0, CUDA.ones(Float32, (M,nvec)),1024, 256)

# Unclear how to solve issue of not having enough threads/blocks when nvec is large. Combining threads and blocks works

## Timing Methods
@btime beta_pdf_gpu($x, $1.0f0, $2.0f0,$1) #35 mus for M =15, nvec = 10000 
@btime beta_pdf_gpu_basic($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)), $1024, $256)  #5.6 s for M =15, nvec = 10000 
@btime beta_pdf_gpu_basic($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)), $1000, $100) #2.3s for M =15, nvec = 10000
@btime beta_pdf_gpu_basic($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)), $1, $1) #360 ms for M =15, nvec = 10000
## Choosing the number of threads and blocks matters a lot in terms of speed, unclear what factors affect the right combination
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$1, $1) # 290 ms for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$16, $16) #1.3 ms for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$128, $16) #210 mus for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$1024, $16) #103 mus for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$1024, $256) #112 mus for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$1024, $1024) #150 mus for M =15, nvec = 10000
@btime beta_pdf_gpu_threads_blocks($x, $1.0f0, $2.0f0, $CUDA.ones(Float32, (M,nvec)),$1024, $2048) #175 mus for M =15, nvec = 10000
# Appears to be that the threads and blocks combo method is fast, but not as fast as the regular way. However, if smart about computing the product, I think it'll be faster 

function beta_pdf_gpu_threads_blocks_prod!(x,a,b,denom,result)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(x)
   #     idx = Int32(ceil(i/M))
        idx = ceil(i/M)
        @inbounds result[idx] *= (x[i]^(a-1.0f0) * (1.0f0-x[i])^(b-1.0f0))/denom
    end
    return nothing
end

function beta_pdf_gpu_threads_blocks_prod(x, a, b, result, threadsPerBlock, blocksPerGrid)
    denom = gamma(a)*gamma(b)/gamma(a+b)
    CUDA.@sync begin
        @cuda threads = threadsPerBlock blocks = blocksPerGrid beta_pdf_gpu_threads_blocks_prod!(x, a, b, denom ,result)
    end
end

beta_pdf_gpu_threads_blocks_prod(x, 1.0f0, 2.0f0, CUDA.ones(Float32, (1,nvec)),1024, 16)
