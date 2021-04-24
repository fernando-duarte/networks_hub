import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA
using SpecialFunctions

@test Threads.nthreads()>1

M= 25 # number of independent uniform random variables
atol=1e-6
rtol=1e-3
nvec=1000000
maxevals=100000000

x = rand(M, nvec)
x32 = convert(Matrix{Float32},x)
cux = CuArray(x)
cux32 = CuArray(x32)

function beta_pdf(x, a, b)
     prod((x.^(a-1.0) .* (ones(M) .- x).^(b-1.0))./(gamma(a)*gamma(b)/gamma(a+b)),dims = 1)
end

function beta_pdf_gpu(x, a, b)
     prod(x.^(a-1.0) .* (CuArray(ones(M)) .- x).^(b-1.0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_ones(x, a, b)
     prod(x.^(a-1.0) .* (ones_mat .- x).^(b-1.0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_ones_32(x, a, b)
     prod(x.^(a-1.0f0) .* (ones_mat .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_ones_32_sync(x, a, b)
    CUDA.@sync prod(x.^(a-1.0f0) .* (ones_mat .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_ones_32_inbounds(x, a, b)
   @inbounds prod(x.^(a-1.0f0) .* (ones_mat .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end


# Benchmarking Different Options
display(@benchmark pdf(Product(Beta.(1.0,2.0*ones(M))),$x)) # running using Beta function from Dsitributions (slowest)
display(@benchmark beta_pdf($x, 1.0,2.0)) #running on cpu using user written Beta function
display(@benchmark beta_pdf($x32, 1.0f0,2.0f0)) #running on cpu using user written Beta function with float 32 (no change)
display(@benchmark beta_pdf_gpu($cux,1.0,2.0)) # running on gpu using user written beta function (faster than cpu)
display(@benchmark beta_pdf_gpu_ones($cux,1.0,2.0)) # running on gpu, user written, preallocated ones (faster than without preallocation)

ones_mat = CuArray(convert(Vector{Float32},ones(M)))
display(@benchmark beta_pdf_gpu_ones_32($cux32,1.0f0,2.0f0)) # running on gpu using user written beta function with float 32 and preallocated ones as float32 (faster)
display(@benchmark beta_pdf_gpu_ones_32_sync($cux32,1.0f0,2.0f0)) # running on gpu using user written beta function with float 32 and preallocated ones as float32 and sync (same as above) 
display(@benchmark beta_pdf_gpu_ones_32_inbounds($cux32,1.0f0,2.0f0)) # running on gpu using user written beta function with float 32 and preallocated ones as float32 and inbounds (same as above) 

