import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA
using SpecialFunctions

@test Threads.nthreads()>1

# User Inputs
M= 10 # number of independent uniform random variables
atol=1e-6
rtol=1e-3
nvec=1000000
maxevals=100000000

# Initializing Matrices
x = CuArray(rand(Float32, (M, nvec)))
ones_mat = CuArray(ones(Float32, M))
x_cpu = rand(Float32, (M, nvec))
ones_mat_cpu = ones(Float32, M)

# Comparing Multiplication Techniques
function parallel_multi(f, x)
   @floop CUDAEx() for i in 1:size(x, 2)
        val = reduce(*,@view(x[:,i]))
        f[i] = val 
    end
    return f
end

result = CUDA.ones(size(x,2),1)
display("Comparing Multiplication Methods")
display(@benchmark parallel_multi(result, $x)) #is 2-5x faster
display(@benchmark reduce(*, $x, dims = 1)) #2-5x slower
display(@benchmark prod($x, dims=1)) #identical to above 

### takeaway: want to use gpu for parallelization of multiplication 

# Comparing Computing PDFs
function beta_pdf_cpu_old(x,a,b)
    pdf(Product(Beta.(a,b*ones(M))),x)
end

function beta_pdf_cpu_new(x, a, b)
     prod(x.^(a-1.0f0) .* (ones_mat_cpu .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu(x, a, b)
     prod(x.^(a-1.0f0) .* (ones_mat .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=1)
end

function beta_pdf_gpu_precompile(x, a, b)
    power1 = a-1.0f0
    power2 = b-1.0f0
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    prod((x.^(power1) .* (ones_mat .- x).^(power2))./denominator,dims=1)
end

function beta_pdf_gpu_precompile_flip(x, a, b)
    power1 = a-1.0f0
    power2 = b-1.0f0
    denominator = exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b))
    ones_mat_flip = ones_mat'
    prod((x.^(power1) .* (ones_mat_flip .- x).^(power2))./denominator,dims=2)
end

function beta_pdf_gpu_precompile_floop(f, x, a, b) # super slow
    power1 = a-1.0f0
    power2 = b-1.0f0
    denominator = (exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b)))[1]
    for j in 1:size(x,2)
        @floop CUDAEx() for w in referenceable(x[:,j])
            @reduce(xprod *= w[]^power1[1] * (1.0f0 - w[])^power2[1]/denominator)
        end
        f[j] = xprod
    end
    return f
end

function beta_pdf_gpu_precompile_floop2(f,x, a, b) # works with CUDAEx() but slightly slower cause of calculations outside of floop loop
    power1 = a-1.0f0
    power2 = b-1.0f0
    denominator = (exp.(CUDA.lgamma.(a)) .* exp.(CUDA.lgamma.(b)) ./ exp.(CUDA.lgamma.(a+b)))[1]
    ones_mat = CuArray(ones(Float32, (M,nvec)))
    val2 = (x.^(power1) .* (ones_mat .- x).^(power2))./denominator
    @floop CUDAEx() for i in 1:size(x,2)
        #f[i] = reduce(*,(x[:,i].^(power1) .* (ones_mat .- x[:,i]).^(power2))./denominator)
        val = reduce(*,@view(val2[:,i]))
        f[i] = val
    end
    return f
end

# Benchmarking PDF Calculations
display("Benchmarking PDF Calculations")
display(@benchmark beta_pdf_cpu_old($x_cpu, 1.0f0, 2.0f0)) #1.6 s M = 25, Nvec = 1000000
display(@benchmark beta_pdf_cpu_new($x_cpu, 1.0f0, 2.0f0)) # 400 ms M = 25, Nvec = 1000000
display(@benchmark beta_pdf_gpu($x, 1.0f0, 2.0f0)) # 33ms M = 25, Nvec = 1000000
display(@benchmark beta_pdf_gpu_precompile($x,1.0f0,2.0f0)) # 28 ms, M = 25 Nvec = 1000000
display(@benchmark beta_pdf_gpu_precompile_flip($x',1.0f0, 2.0f0)) # 28ms M = 25, Nvec = 1000000

result = CUDA.ones(Float32, (size(x,2),1))
display(@benchmark beta_pdf_gpu_precompile_floop(result, $x,1.0f0, 2.0f0)) # 180 seconds for M = 25, Nvec = 1000000
display(@benchmark beta_pdf_gpu_precompile_floop2(result, $x, 1.0f0, 2.0f0)) # 90 ms for M = 25, Nvec  = 1000000

# Benchmarking PDF Calculations
display("Comparing PDF Calculations")
display(@benchmark beta_pdf_cpu($x_cpu, 1.0f0, 2.0f0))
display(@benchmark beta_pdf_gpu($x, 1.0f0, 2.0f0)) #10-15x faster than cpu
display(@benchmark beta_pdf_gpu_precompile($x,1.0f0,2.0f0)) #same as above 
result = CUDA.ones(Float32, (size(x,2),1))
display(@benchmark beta_pdf_gpu_precompile_floop($result, $x,1.0f0, 2.0f0)) #a bit slower than aboveh slower

### Takeaway: Want to parallelize the pdf computation using floop with GPU, but have to figure out how to also do calculations with it to get gain from multiplication
