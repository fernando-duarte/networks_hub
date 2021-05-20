import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test, CUDA
using FLoops, FoldsCUDA
using SpecialFunctions
using Suppressor

@test Threads.nthreads()>1

# User Inputs
M= 5 # number of independent uniform random variables
atol=1e-10
rtol=1e-10
nvec=100000
maxevals=1000000

# Initializing Matrices
ones_mat = CuArray(ones(Float32, M))
result = CUDA.ones(Float32, (nvec,1))
x_cpu = rand(Float64, (M, nvec))
x = convert(CuArray{Float32},x_cpu)
x_1d_cpu = rand(Float64, M)
x_1d_gpu = convert(CuArray{Float32},x_1d_cpu)

#Initializing Timers
time_convert = 0
time_calc_cpu = 0
time_calc_gpu = 0

# Initializing pdf functions
function beta_pdf_gpu(x::CuArray,a,b,dim)
    prod(x.^(a-1.0f0) .* (1.0f0 .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=dim)
end

function beta_pdf_cpu(x,a,b,dim)
    prod(x.^(a-1.0f0) .* (1.0f0 .- x).^(b-1.0f0)./(gamma(a)*gamma(b)/gamma(a+b)),dims=dim)
end

# Timing pdf functions 
display(@benchmark convert(CuArray{Float32}, $x_cpu)) #500 mus for M = 5 and Nvec = 100000

display(@benchmark beta_pdf_gpu($x, $1.0f0, $2.0f0, $1))  #2.3 ms for M = 5 and Nvec = 100000, 170 mus on p3
display(@benchmark beta_pdf_cpu($x_cpu, $1.0f0, $2.0f0, $1)) #7.3 ms for M = 5 and Nvec = 100000

display(@benchmark beta_pdf_gpu($x[:,1], $1.0f0, $2.0f0, $1)) #50 mus for M = 5 and Nvec = 100000
display(@benchmark beta_pdf_cpu($x_cpu[:,1], $1.0f0, $2.0f0, $1)) #400 ns for M = 5 and Nvec = 100000

display(@benchmark beta_pdf_gpu($x_1d_gpu, $1.0f0, $2.0f0, $1)) #40 mus for M = 5 and Nvec = 100000
display(@benchmark beta_pdf_cpu($x_1d_cpu, $1.0f0, $2.0f0, $1)) #300 ns for M = 5 and Nvec = 100000

@test beta_pdf_gpu(x, 1.0f0, 2.0f0, 1) ≈ convert(CuArray{Float32},beta_pdf_cpu(x_cpu, 1.0f0, 2.0f0, 1))
@test beta_pdf_gpu(x[:,1], 1.0f0, 2.0f0, 1)[1] ≈ beta_pdf_cpu(x_cpu[:,1], 1.0f0, 2.0f0, 1)[1]
@test beta_pdf_gpu(x_1d_gpu, 1.0f0, 2.0f0, 1)[1] ≈ beta_pdf_cpu(x_1d_cpu, 1.0f0, 2.0f0, 1)[1]

## Takeaway: For large nvec, gpu calculations are faster even when factoring in the time to convert to gpu array (2.3ms + .4ms < 7.3ms). However, for small nvec (nvec == 1) cpu calculation is faster. At around nvec = 500 the gpu and cpu are equalivalently fast. The second example shows how we do it for the col for loop, and that the cpu is much faster. We can't compute column by column with gpu, have to compute all at once for gpu to have speed advantage

## Integration Functions
function cpu_2_gpu(x,f)
    global time_convert += @elapsed x = convert(CuArray{Float32},x)
    global time_calc_gpu += @elapsed f[1] = beta_pdf_gpu(x, 1.0f0, 2.0f0, 1)[1]
end

function cpu_2_gpu_col(x, f)
   global time_convert += @elapsed  x = convert(CuArray{Float32},x)
    Threads.@threads for i in 1:size(x,2)
       global time_calc_gpu += @elapsed f[i] = beta_pdf_gpu(x[:,i], 1.0f0, 2.0f0,1)[1]
    end
end


function cpu_2_gpu_col2(x, f)
   global time_convert += @elapsed x = convert(CuArray{Float32},x)
    global time_calc_gpu += @elapsed val = beta_pdf_gpu(x, 1.0f0, 2.0f0,1)
    Threads.@threads for i in 1:size(x,2)
       @elapsed f[i] = val[i]
    end
end

## Cpu Comparison

function int_cpu(x, f)
  global time_calc_cpu += @elapsed f[1] = beta_pdf_cpu(x,1.0,2.0,1)[1]
end

function int_thread_col_cpu(x, f)
    @floop for i in 1:size(x,2)
        global time_calc_cpu += @elapsed f[i] = beta_pdf_cpu(x[:,i],1.0,2.0,1)[1]
    end
end

# Integration
result_gpu = cuhre(cpu_2_gpu, M, 1, atol=atol, rtol=rtol, nvec = 1)
result_cpu = cuhre(int_cpu, M, 1, atol=atol, rtol=rtol,nvec = 1)
display(time_calc_cpu)
display(time_calc_gpu)
display(time_convert)

@test abs(result_gpu[1][1] - result_cpu[1][1]) < 1e-5 #errors but essentially equal

#Takeaway: Converting cost is too high for 1d case. In addition, the calculations are longer for gpu than cpu because nvec = 1 in this case

time_calc_cpu = 0 
time_calc_gpu = 0
time_convert = 0

@suppress begin #blocking warnings
#    suave(cpu_2_gpu_col, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals)
    result_gpu = suave(cpu_2_gpu_col2, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals)
    result_cpu = suave(int_thread_col_cpu, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals)
end

display(time_calc_cpu)
display(time_calc_gpu)
display(time_convert)

@test abs(result_gpu[1][1] - result_cpu[1][1]) < 1e-5 #errors but approx equal

# Takeaways: with larger nvec conversion time + gpu time is below cpu calc, 


## Integration Functions
function cpu_2_gpu(x,f)
    x = convert(CuArray{Float32},x)
    f[1] = beta_pdf_gpu(x, 1.0f0, 2.0f0, 1)[1]
end

function cpu_2_gpu_col(x, f)
   x = convert(CuArray{Float32},x)
    Threads.@threads for i in 1:size(x,2)
       f[i] = beta_pdf_gpu(x[:,i], 1.0f0, 2.0f0,1)[1]
    end
end


function cpu_2_gpu_col2(x, f)
   x = convert(CuArray{Float32},x)
    val = beta_pdf_gpu(x, 1.0f0, 2.0f0,1)
    Threads.@threads for i in 1:size(x,2)
       @elapsed f[i] = val[i]
    end
end

## Cpu Comparison

function int_cpu(x, f)
    f[1] = beta_pdf_cpu(x,1.0,2.0,1)[1]
end

function int_thread_col_cpu(x, f)
    @floop for i in 1:size(x,2)
        f[i] = beta_pdf_cpu(x[:,i],1.0,2.0,1)[1]
    end
end

# Integration
display(@benchmark cuhre($cpu_2_gpu, $M, $1, atol=$atol, rtol=$rtol)) # 100s for M = 5, 80s on p3
@suppress begin
    display(@benchmark suave($cpu_2_gpu_col, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 169 s for M = 5, 141s for p3
    display(@benchmark suave($cpu_2_gpu_col2, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 42 s for M = 5, 26s for p3
end

display(@benchmark cuhre($int_cpu, $M, $1, atol=$atol, rtol=$rtol)) # 350 ms for M = 5
display(@benchmark suave($int_thread_col_cpu, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #800ms for M =5  

@suppress begin
    result_gpu = suave(cpu_2_gpu_col2, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    result_cpu = suave(int_thread_col_cpu, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    @test result_gpu[1][1] - result_cpu[1][1] < 1e-5
end


## Takeaway CPU looks much faster unless we can increase the size of nvec drastically. 

# Trying with larger M 
M = 15
maxevals = 10000000
@suppress begin
    display(@benchmark suave($cpu_2_gpu_col2, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) # 250s for p3
    display(@benchmark suave($int_thread_col_cpu, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals)) #19s for p3
    result_gpu = suave(cpu_2_gpu_col2, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    result_cpu = suave(int_thread_col_cpu, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    @test result_gpu[1][1] - result_cpu[1][1] < 1e-2 # but don't equal 1...?
end

M = 25
@suppress begin
    display(@benchmark suave($cpu_2_gpu_col2, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals))
    display(@benchmark suave($int_thread_col_cpu, $M, $1, atol=$atol, rtol=$rtol, nvec = $nvec, maxevals = $maxevals))
    result_gpu = suave(cpu_2_gpu_col2, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    result_cpu = suave(int_thread_col_cpu, M, 1, atol=atol, rtol=rtol, nvec = nvec, maxevals = maxevals) 
    @test result_gpu[1][1] - result_cpu[1][1] < 1e-2 # but don't equal 1...?
end
