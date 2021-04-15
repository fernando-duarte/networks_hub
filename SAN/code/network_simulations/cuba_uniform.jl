import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()
using Cuba, Distributions
using BenchmarkTools, Test
using CSV, DataFrames

@test Threads.nthreads()>1

M= 5 # number of independent beta random variables
atol=1e-6
rtol=1e-3
nvec = 100000

# integrate the pdf of the joint distribution -- should always equal 1
function int(x, f)
   f[1] = pdf(Product(Uniform.(zeros(M),ones(M))),x)
end

# multithread
function int_thread_col(x, f)
    Threads.@threads for i in 1:size(x,2)
      f[i] = pdf(Product(Uniform.(zeros(M),ones(M))),@view(x[:,i]))
    end
end

# multithread and loop to create product distribution
function int_thread_el(x,f)
   f[1,:] .= 1.0
   Threads.@threads for j in 1:size(x,2)
       for i in 1:size(x, 1)
           f[1, j] *= pdf(Uniform(0.0,1.0),@view(x[i,j]))
       end
   end
end

# Cuhre
time_int = @benchmark cuhre($int, $M, 1, atol=$atol, rtol=$rtol)
result_int = cuhre(int, M, 1, atol=atol, rtol=rtol)[1]
time_col = @benchmark cuhre($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec)
result_col = cuhre(int_thread_col, M, 1, atol=atol, rtol=rtol, nvec=nvec)[1]
time_el = @benchmark cuhre($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec)
result_el = cuhre(int_thread_el, M, 1, atol=atol, rtol=rtol, nvec=nvec)[1]

data_out = DataFrame(M = M, alg = "cuhre", atol = atol, rtol = rtol, nvec = 100000, max_evals = 0,
    mean_int = result_int, mean_col = result_col, mean_el = result_el, 
    time_int_sec = median(time_int).time/1e9, time_col_sec = median(time_col).time/1e9, time_el_sec = median(time_el).time/1e9,
    mem_int_gb = time_int.memory/1000/2^20, mem_col_gb = time_col.memory/1000/2^20, mem_el_gb = time_el.memory/1000/2^20)
CSV.write("data_out_$(M)_cuhre.csv",  data_out)

# Divonne
maxevals=1174707384

time_int = @benchmark divonne($int, $M, 1, atol=$atol, rtol=$rtol,maxevals=$maxevals)
result_int = divonne(int, M, 1, atol=atol, rtol=rtol,maxevals=maxevals)[1]
time_col = @benchmark divonne($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec,maxevals=$maxevals)
result_col = divonne(int_thread_col, M, 1, atol=atol, rtol=rtol, nvec=nvec,maxevals=maxevals)[1]
time_el = @benchmark divonne($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec,maxevals=$maxevals)
result_el = divonne(int_thread_el, M, 1, atol=atol, rtol=rtol, nvec=nvec,maxevals=maxevals)[1]

data_out = DataFrame(M = M, alg = "divonne", atol = atol, rtol = rtol, nvec = 100000, max_evals = maxevals,
    mean_int = result_int, mean_col = result_col, mean_el = result_el, 
    time_int_sec = median(time_int).time/1e9, time_col_sec = median(time_col).time/1e9, time_el_sec = median(time_el).time/1e9,
    mem_int_gb = time_int.memory/1000/2^20, mem_col_gb = time_col.memory/1000/2^20, mem_el_gb = time_el.memory/1000/2^20)
CSV.write("data_out_$(M)_divonne.csv",  data_out)


# monte carlo suave
nmin=2
nnew=80000
flatness=150

time_int = @benchmark suave($int, $M, 1, atol=$atol, rtol=$rtol,maxevals=$maxevals,nnew=$nnew,nmin=$nmin,flatness=$flatness)
result_int = suave(int, M, 1, atol=atol, rtol=rtol,maxevals=maxevals,nnew=nnew,nmin=nmin,flatness=flatness)[1]
time_col = @benchmark suave($int_thread_col, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec,maxevals=$maxevals,nnew=$nnew,nmin=$nmin,flatness=$flatness)
result_col = suave(int_thread_col, M, 1, atol=atol, rtol=rtol, nvec=nvec,maxevals=maxevals,nnew=nnew,nmin=nmin,flatness=flatness)[1]
time_el = @benchmark suave($int_thread_el, $M, 1, atol=$atol, rtol=$rtol,nvec=$nvec,maxevals=$maxevals,nnew=$nnew,nmin=$nmin,flatness=$flatness)
result_el = suave(int_thread_el, M, 1, atol=atol, rtol=rtol, nvec=nvec,maxevals=maxevals,nnew=nnew,nmin=nmin,flatness=flatness)[1]

data_out = DataFrame(M = M, alg = "suave", atol = atol, rtol = rtol, nvec = 100000, max_evals = maxevals,
    mean_int = result_int, mean_col = result_col, mean_el = result_el, 
    time_int_sec = median(time_int).time/1e9, time_col_sec = median(time_col).time/1e9, time_el_sec = median(time_el).time/1e9,
    mem_int_gb = time_int.memory/1000/2^20, mem_col_gb = time_col.memory/1000/2^20, mem_el_gb = time_el.memory/1000/2^20)
CSV.write("data_out_$(M)_suave.csv",  data_out)

# ### Default values of parameters
# # Common arguments.
# const NVEC      = 1
# const RTOL      = 1e-4
# const ATOL      = 1e-12
# const FLAGS     = 0
# const SEED      = 0
# const MINEVALS  = 0
# const MAXEVALS  = 1000000
# const STATEFILE = ""
# const SPIN      = C_NULL

# # Vegas-specific arguments.
# const NSTART    = 1000
# const NINCREASE = 500
# const NBATCH    = 1000
# const GRIDNO    = 0

# # Suave-specific arguments.
# const NNEW     = 1000
# const NMIN     = 2
# const FLATNESS = 25.0

# # Divonne-specific arguments.
# const KEY1         = 47
# const KEY2         = 1
# const KEY3         = 1
# const MAXPASS      = 5
# const BORDER       = 0.0
# const MAXCHISQ     = 10.0
# const MINDEVIATION = 0.25
# const NGIVEN       = 0
# const LDXGIVEN     = 0
# const XGIVEN       = 0
# const NEXTRA       = 0
# const PEAKFINDER   = C_NULL

# # Cuhre-specific argument.
# const KEY = 0
