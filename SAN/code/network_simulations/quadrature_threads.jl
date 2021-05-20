import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using Cuba, Distributions
using BenchmarkTools, Test

# Inputs
M= 15 # number of independent beta random variables
atol=1e-6
rtol=1e-3

global num_evals = 1 # shows x if 1, doesn't otherwise

# integrate the pdf of the joint distribution -- should always equal 1
function int(x, f)
    if num_evals == 1 
        @show(size(x))
    end
    global num_evals += 1
   f[1] = pdf(Product(Beta.(1.0,2.0*ones(M))),x)
end

# multithread
function int_thread_col(x, f)
     if num_evals == 1 
        @show(size(x))
    end
    global num_evals += 1
    Threads.@threads for i in 1:size(x,2)
      f[i] = pdf(Product(Beta.(1.0,2.0*ones(M))),@view(x[:,i]))
    end
end

# multithread and loop to create product distribution
function int_thread_el(x,f)
     if num_evals == 1 
        @show(size(x))
    end
    global num_evals += 1
   f[1,:] .= 1.0
   Threads.@threads for j in 1:size(x,2)
       for i in 1:size(x, 1)
           f[1, j] *= pdf(Beta(1.0,2.0),@view(x[i,j]))
       end
   end
end

# we get the right answer
@show result, err = cuhre(int, M, 1, atol=atol, rtol=rtol)
@show result, err = cuhre(int_thread_col, M, 1, atol=atol, rtol=rtol,nvec=100)
@show result, err = cuhre(int_thread_el, M, 1, atol=atol, rtol=rtol,nvec=100)

@btime cuhre($(int), $M, 1, atol=$atol, rtol=$rtol) # slow

println("MultiThreading with Cuhre")
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(16)) #size is 16 for M = 5, size is 16 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100)) #size is 100 for M = 5, size is 100 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000)) #size is 273 for M = 5, size is 1000 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000)) #size is 273 for M = 5, size is 10000 for M = 15, 
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100000)) #size is 273 for M = 5, size is 37789 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000000)) #size is 273 for M = 5, size is 37789 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000000)) #size is 273 for M = 5,size is 37789 for M = 15,

println("MultiThreading with Suave")
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(16)) #size is 16 for M = 5, size is 16 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100)) #size is 100 for M = 5, size is 100 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000)) #size is 1000 for M = 5, size is 1000 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000)) #size is 1000 for M = 5, size is 1000 for M = 15, 
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100000)) #size is 1000 for M = 5, size is 1000 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000000)) #size is 1000 for M = 5, size is 1000 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000000)) #size is 1000 for M = 5,size is 1000 for M = 15,


println("MultiThreading with Divonne")
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(16)) #size is 24 for M = 5, size is 16 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100)) #size is 24 for M = 5, size is 24 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000)) #size is 24 for M = 5, size is 24 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000)) #size is 24 for M = 5, size is 24 for M = 15, 
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100000)) #size is 24 for M = 5, size is 24 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000000)) #size is 24 for M = 5, size is 24 for M = 15,
global num_evals = 1 #resetting so it shows size(x)
@btime divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000000)) #size is 24 for M = 5,size is 24 for M = 15,

println("multithread and create product in loop")
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(16)) #size is 16 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100)) #size is 100 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000)) #size is 273 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000)) #size is 273 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(100000)) #size is 273 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(1000000)) #size is 273 for M = 5
global num_evals = 1 #resetting so it shows size(x)
@btime cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$(10000000)) #size is 273 for M = 5
global num_evals = 1 #resetting so it shows size(x)

## Int_thread_el and Int_thread_col have same size(x) for given M and nvec
## Takeaway speed increases with nvec for cuhre because the size of the x matrix grows from 16 to 30k for M = 15. 
## Takeaway, Speed seems to increase from nvec = 16 to nvec = 100 to nvec = 1000 for Suave and Divonne, but after that remains constant as size of X maxes out at like 250. 
