import Pkg; Pkg.add("Cbc")
import Pkg; Pkg.add("Juniper")
import Pkg; Pkg.add("JuMP")
import Pkg; Pkg.add("InfiniteOpt")
import Pkg; Pkg.add("Ipopt")
import Pkg; Pkg.add("SCIP")
import Pkg; Pkg.add("XLSX")
import Pkg; Pkg.add("DataFrames")
import Pkg; Pkg.add("Random")

using InfiniteOpt, JuMP, Distributions, LinearAlgebra, Random
using XLSX, DataFrames
using Juniper, Ipopt, Cbc, SCIP

N = 5 # keep largest `N' nodes by assets

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units

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
M = 1000 #large M

# initial guess
rng =  Random.seed!(123)
A0 = zeros(N,N)
c0 =  data_nm.c 
b0 = data_nm.b
#α0 = 1.0*ones(N)
#β0= 50.0*ones(N) 

## Plugging in answers from benchmark_comp
α0 = [1.9769572880132678, 2.03241259121659, 1.9932613374527037, 1.87706227793852, 1.6069831835195945] 
β0 = [52.560963982981754, 40.34375025060937, 47.902760862643134, 78.19354871489449, 104.17269222248204] 
c0 = [1.451636125, 1.38565475, 1.516806125, 1.158976375, 1.0802461084137611] 
b0 = [1.2761951138380403, 1.7885790392572256, 1.6398795667566521, 1.206356349290088, 0.039981914271755244] 
A0 = [0.0 0.14133583887525003 0.15110396086528224 0.07168597756011581 1.2421834356036085e-45; 1.6219083799719297e-46 1.2005708277365943e-46 0.001113034911629969 0.0021575868855880937 6.201883233297931e-47; 0.0 0.0 1.155259280558607e-45 0.0017754459191141763 0.000340912355544804; 0.0 0.0 0.0 2.039381281484266e-46 0.0008006562534731075; 0.7006089517923049 0.26066963219151473 6.353882944741388e-46 0.0 0.0] 

# Initialize the  model 

m = InfiniteModel(SCIP.Optimizer)

# Define the random parameters
@infinite_parameter(m, x[1:N] in Beta(1.0, 50.0), num_supports = 1000) 

# Define the variables
@infinite_variable(m, p[1:N](x)) # clearing vector. 
@infinite_variable(m, z1[1:N](x))

# Define A as a matrix of variables: 
@hold_variable(m, c[i=1:N]  >= 0) #outside assets for node i, 
@hold_variable(m, b[i=1:N] >= 0) #outside liabilities for node i
@hold_variable(m, A[i=1:N, j=1:N] >= 0) #net payments scaled by total liabilites for firm i between nodes i and j
@hold_variable(m, α[i=1:N] >= 0) # parameter that determines the shocks
@hold_variable(m, β[i=1:N] >= 0) # parameter that determines shocks 
@constraint(m, w .<= c .<= data.assets) # outside assets greater than networth and lower than total assets
@constraint(m, 0 .<= b .<= data.p_bar) #outside liabilities greater than 0 and less than total liabilities
@constraint(m, 0.8 .<= α .<= 3.0) 
@constraint(m, 1.0 .<= β .<= 150.0)
@constraint(m, -sum(A,dims=2).*data.p_bar .+ data.p_bar .==  b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d


# Constraints for A
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

# Fixing variables
@constraint(m, α .== α0)
@constraint(m, β .== β0)
@constraint(m, b .== b0)
@constraint(m, c .== c0)
@constraint(m, A .== A0)

# Define the objective
@objective(m, Min, expect(sum(-p),x))

# Define the constraints for 
# p == min(p_bar, max(0, (1+g0).* A'*p + c - x.*c - g0*p_bar)
@infinite_variable(m, ind[1:N](x), Bin)
@infinite_variable(m, ind2[1:N](x), Bin)
@constraint(m, p .<= p_bar)
@constraint(m, p .<= (1+g0) .* A * p - x.*c .- g0.*p_bar)
@constraint(m, p .>= p_bar - M .* ind)
@constraint(m, p .>= (1+g0) .* A * p - x.*c .- g0.*p_bar - M .* (1 .- ind))

@constraint(m, z1 .>= zeros(N))
@constraint(m, z1 .>=  (1+g0) .* A * p - x.*c .- g0.*p_bar)
@constraint(m, z1 .<= zeros(N) + M .* ind2)
@constraint(m, z1 .<= (1+g0) .* A * p - x.*c .- g0.*p_bar + M .* (1 .- ind2))

# Solve the model and get the results
optimize!(m)
opt_objective = objective_value(m)
opt_A = value.(A)