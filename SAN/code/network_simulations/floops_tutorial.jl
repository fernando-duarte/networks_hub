using FLoops, BenchmarkTools
using FoldsCUDA, CUDA
using GPUArrays: @allowscalar

function demo(executor)
    s = 1.0
           @floop executor for x in 1:6
               @reduce(s *= x)
           end
           return s
       end;

@benchmark demo(SequentialEx(simd = Val(true))) #sequentialEX seems way faster than other options surprisingly 
@benchmark demo(ThreadedEx(basesize = 2))
@benchmark demo(DistributedEx(threads_basesize = 2))

# Compute PI 
xs = rand(10000)
ys = rand(10000)
p = count(xs.^2 .+ ys.^2 .< 1) / length(xs)
4 * p

using Random123

rng_a = Philox2x(0)
rng_b = Philox2x(0)
set_counter!(rng_b, typemax(UInt64))
rand(rng_b, UInt64, 2)
@assert rng_a == rng_b

function counters(n)
    stride = typemax(UInt64) ÷ n
    return UInt64(0):stride:typemax(UInt64)-stride
end

function monte_carlo_pi(n, m = 10_000, ex = ThreadedEx())
    @floop ex for ctr in counters(n)
        rng = set_counter!(Philox2x(0), ctr)
        nhits = 0
        for _ in 1:m
            x = rand(rng)
            y = rand(rng)
            nhits += x^2 + y^2 < 1
        end
        @reduce(tot = 0 + nhits)
    end
    return 4 * tot / (n * m)
end

function monte_carlo_pi_gpu(n, m = 10_000, ex = has_cuda_gpu() ? CUDAEx() : ThreadedEx())
    @floop ex for ctr in counters(n)
        rng = set_counter!(Philox2x(0), ctr)
        nhits = 0
        for _ in 1:m
            x = rand(rng)
            y = rand(rng)
            nhits += x^2 + y^2 < 1
        end
        @reduce(tot = 0 + nhits)
    end
    return 4 * tot / (n * m)
end

@benchmark monte_carlo_pi(2^6)
@benchmark monte_carlo_pi_gpu(2^6) #out of gpu memory error for larger N 

## Find max from an array - so much faster on GPU
xs = CUDA.rand(10^4);

function findmax_gpu()
 @floop CUDAEx() for (x, i) in zip(xs, eachindex(xs))
       @reduce() do (imax = -1; i), (xmax = -Inf32; x)
          if xmax < x
             xmax = x
             imax = i
          end
      end
end
end
    
@benchmark findmax_gpu()

xs = rand(10^4)
function findmax_cpu()
 @floop  for (x, i) in zip(xs, eachindex(xs))
      @reduce() do (imax = -1; i), (xmax = -Inf32; x)
          if xmax < x
             xmax = x
             imax = i
          end
     end
 end
end

@benchmark findmax_cpu()
