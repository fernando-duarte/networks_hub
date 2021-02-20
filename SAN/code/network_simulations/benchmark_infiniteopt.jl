using Pkg
Pkg.add(Pkg.PackageSpec(;name="InfiniteOpt", version="0.3.2"))

using LinearAlgebra, DataFrames, XLSX, JuMP, Ipopt, Distributions, Random
using InfiniteOpt, Test

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

# initial guess
rng =  Random.seed!(123)
A0 = [0.0 0.14133583887525153 0.1511039608653196 0.07168597756012116 0.0; 4.397790928630315e-46 3.784129362573653e-46 0.0011130349115881636 0.0021575868856038793 0.0; 0.0 2.6904824890119278e-46 0.0 0.0017754459190904477 0.0003409123555305989; 0.0 8.151907446024649e-47 4.0774174184842646e-46 1.3738471521896145e-46 0.0008006562534392179; 0.7006089517923049 0.26066963219151196 1.0717105137003107e-47 0.0 0.0] #zeros(N,N)
c0 = [1.451636125, 1.38565475, 1.516806125, 1.158976375, 1.0802461084138253] # data_nm.c 
b0 = [1.2761951138379517, 1.7885790392572722, 1.6398795667567143, 1.206356349290129, 0.03998191427175805] #data_nm.b
#α0 = 1.0*ones(N)
#β0= 50.0*ones(N) 
α0 = [1.9769572880132678, 2.0324125912165893, 1.993261337452703, 1.8770622779385235, 1.6069831835194497] 
β0 = [52.560963982981754, 40.34375025060936, 47.90276086264312, 78.19354871489458, 104.1726922224828] 
α0 = α0[1:N]
β0 = β0[1:N]

# Setting up model 
m = InfiniteModel(Ipopt.Optimizer);

# Parameter everything depends upon (shocks to each node)
@infinite_parameter(m, x[i=1:N] in Beta(α0[i],β0[i]), num_supports = 10000) #parameter of shocks

# setting up variables, starting positino, and constraints
@infinite_variable(m, p[i=1:N](x), start = 0) # clearing vector. 
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

# Fixing variables
@constraint(m, α .== α0)
@constraint(m, β .== β0)
@constraint(m, b .== b0)
@constraint(m, c .== c0)

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

# Clearing Vector Constraints
#@constraint(m, p == min(p_bar, max(zeros(N), (1+g0) .* transpose(A)*p + c - x .* c .- g0.*p_bar))) 
@constraint(m, p .== (1+g0) .* transpose(A)*p + c - x .* c .- g0.*p_bar)

# Distributional constraints
#delta = Pr(x*c >= w)
#@infinite_variable(m, y[i=1:N](x), Bin)
#@constraint(m, w .- c .* x .<= y .* M)
#for i = 1:N
#    @constraint(m, expect(1 .- y[i], x[i]) == delta[i])
#end

# Setting up objective function
@objective(m, Min, expect(sum(c.*x + p_bar - p),x)) # min/max E[(ci * xi + p_bar i - pi(x)))]

# Training Model
@time optimize!(m);

## Checking Solution 
# variable values after training
csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)
psol0 = JuMP.value.(p)
objective = objective_value(m)
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

# Displaying output in batch solution
print("c solution = $csol0 \n")
print("c initial = $c0 \n")
print("b solution = $bsol0 \n")
print("b initial = $b0 \n")
print("A solution = $Asol0 \n")
print("A initial = $A0 \n")
print("α solution = $αsol0 \n")
print("α initial = $α0 \n")
print("β solution = $βsol0 \n")
print("β initial = $β0 \n")
print("Final value = $objective \n")