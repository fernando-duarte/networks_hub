using Pkg
Pkg.build("HDF5")
Pkg.add("JLD")
Pkg.build("SpecialFunctions")
Pkg.add("SpecialFunctions")
Pkg.build("FFTW")
Pkg.add("FFTW")
Pkg.add("PATHSolver")
Pkg.build("SLEEFPirates")
Pkg.add("SLEEFPirates")
Pkg.add("IterTools")
Pkg.add("Tracker")
Pkg.add("InfiniteOpt")
Pkg.build("DiffEqBase")
Pkg.add("DiffEqBase")
Pkg.add("DistributionsAD")

# restart kernel after running the lines above
using Pkg
Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("Complementarity")
Pkg.add("Zygote")
Pkg.add(Pkg.PackageSpec(url="https://github.com/SciML/Quadrature.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaDiff/ForwardDiff.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaDiff/FiniteDiff.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/giordano/Cuba.jl", rev="master"))
Pkg.add("Cubature")
Pkg.add("HCubature")
Pkg.add("PlotlyBase");Pkg.add("PlotlyJS"); Pkg.add("ORCA")
Pkg.add("NLopt")
 
using LinearAlgebra, DataFrames, XLSX, Missings, JuMP, Ipopt, Random, Complementarity, Test, Distributions, DistributionsAD, SpecialFunctions, NLsolve
using Quadrature, ForwardDiff, FiniteDiff, Zygote, Cuba, Cubature, HCubature
using Plots, InfiniteOpt
plotly(ticks=:native)  

N = 10 # keep largest `N' nodes by assets

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units
# data.b[:] .= missing
# data.c[:] .= missing

col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing
data_nm = coalesce.(data, data.assets/1.5) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))
dropmissing(data, [:delta, :delta_alt, :w, :assets, :p_bar]) # remove type missing

names(data) # column names
describe(data)
show(data, allcols = true)

p_bar = data.p_bar
assets = data.assets
delta = data.delta
w= data.w
g0 = 0.05   # bankruptcy cost

# initial guess
rng =  Random.seed!(123)
A0 = zeros(N,N)
c0 =  data_nm.c 
b0 = data_nm.b
α0 = 1.0*ones(N)
β0= 50.0*ones(N) 

# Setting up model 
m = InfiniteModel(Ipopt.Optimizer);

# Parameter everything depends upon (shocks to each node)
@infinite_parameter(m, x[i=1:N] in Beta(α0[i],β0[i])) #parameter of shocks

# setting up variables, starting positino, and constraints
@infinite_variable(m, p[i=1:N](x), start = 0) # clearing vector. 
@hold_variable(m, c[i=1:N]) #outside assets for node i, 
@hold_variable(m, b[i=1:N]) #outside liabilities for node i
@hold_variable(m, A[i=1:N, j=1:N]) #net payments scaled by total liabilites for firm i between nodes i and j
@hold_variable(m, α[i=1:N]) # parameter that determines the shocks
@hold_variable(m, β[i=1:N]) # parameter that determines shocks 

for i = 1:N
    @constraint(m, w[i] <= c[i] <= data.assets[i]) # outside assets greater than networth and lower than total assets
    @constraint(m, 0 <= b[i] <= data.p_bar[i]) #outside liabilities greater than 0 and less than total liabilities
    @constraint(m, 0.8 <= α[i] <= 3.0)
    @constraint(m, 1.0 <= β[i] <= 150.0)
end


@constraint(m, -sum(A,dims=2).*data.p_bar .+ data.p_bar .==  b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d

for i = 1:N
    for j = 1:N
        @constraint(m, 0<= A[i, j] <= 1)
    end
end

for i = 1:N
    @constraint(m, A[i,i]==0)
end

for i=1:N
    j=1
    while j < i
        @constraint(m, 0 == A[i,j] * A[j,i])
        j += 1
    end
end

# Clearing Vector Constraints
#@constraint(m, p == min(p_bar, max(zeros(N), transpose(A)*p + c - x))) 

@constraint(m, p .== transpose(A)*p + c - x)

## These constraitns were working but then broke
#for i = 1:N
#    @BDconstraint(m, ((transpose(A)*p+c-x)[i] in [-1000000.0,-1e-16]), p[i] == 0) # transpose(A) * p + c -x <= 0 
#    @BDconstraint(m, ((transpose(A)*p+c-x)[i] in [0.0, p_bar[i]]), p[i] == (transpose(A)*p+c-x)[i])
#    @BDconstraint(m, ((transpose(A)*p+c-x)[i] in [p_bar[i]+1e-9, 1000000.0]), p[i] == p_bar[i])
#end

# Distributional constraints
#delta = Pr(x*c >= w)
M = 10000
@infinite_variable(m, y[i=1:N](x), Bin)
@constraint(m, w .- c .* x .<= y .* M)
for i = 1:N
    @constraint(m, expect(1 .- y[i], x[i]) == delta[i])
end


# Setting up objective function
@objective(m, Max, expect(sum(c.*x + p_bar - p),x)) # min/max E[(ci * xi + p_bar i - pi(x)))]

# Training Model
@time JuMP.optimize!(m);

## Checking Solution 
# variable values after training
csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)
tol = 1e-5

#tests 
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(-tol .<=Asol0.<=1+tol)
    @test all(-tol .<=bsol0.<=data.p_bar)
    @test all(-tol .<=csol0.<=data.assets)   
end
