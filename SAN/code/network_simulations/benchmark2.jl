using COSMO
using UnicodePlots
using Convex, SCS
using SparseArrays 
using Distributions
using ProxSDP, CSDP

using LinearAlgebra, DataFrames, XLSX, Missings, JuMP, Ipopt, Random, Complementarity, LinearAlgebra, Test, Expectations

N = 10 # keep largest `N' nodes by assets

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... )
unique!(data) # delete duplicate rows, use `nonunique(data)' to see if there are any duplicates
data = data[isequal.(data.qt_dt,192), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
#N = size(data,1) # number of nodes
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units
data.b[:] .= missing
data.c[:] .= missing

col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing
data_nm = coalesce.(data, data.w .+ 0.01) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))
dropmissing(data, [:delta, :delta_alt, :w, :assets, :p_bar]) # remove type missing

names(data) # column names
describe(data)
show(data, true)

p_bar = data_nm.p_bar
assets = data_nm.assets
delta = data_nm.delta
w= data_nm.w

g0 = 0.0   # bankruptcy cost
rng =  Random.seed!(123)
α=1.0f0 ;β=2.0f0 ;
D = 1 # number of draws
x0 = rand(N,D)
A0 = rand(N,N)

# set up optimization
# m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps",max_iter=10000,hessian_approximation="limited-memory"))
m = Model(with_optimizer(Ipopt.Optimizer,max_iter=10000,print_level=5))

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = data_nm.c[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = data_nm.b[i])   

@constraint(m, sum(A,dims=2).*data.p_bar .== data.p_bar .- b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d

# liabilities are net liabilities: A[i,i]=0 and A[i,j]A[j,i]=0
@constraint(m, [i = 1:N], A[i,i]==0)
for i=1:N
    j=1
    while j < i
        @complements(m, 0 <= A[i,j],  A[j,i] >= 0)
        j += 1
    end
end

# register max and min non-linear functions
maxfun(n1, n2) = max(n1, n2)
minfun(n1, n2) = min(n1, n2)
JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# clearing vector
myexpr = []
for j=1:D
    push!(myexpr, (1+g0).*(A'*p[:,j] .+ c .- x0[:,j].*c ) .- g0.*data.p_bar)
end
    
@variable(m, aux[i=1:N,j=1:D])
for i = 1:N
    for j=1:D
        @constraint(m, aux[i,j] .== myexpr[j][i] )
    end
end

for i = 1:N
    for j=1:D
        @NLconstraint(m,  minfun(data.p_bar[i], maxfun(aux[i,j],0)) == p[i,j] ) 
    end
end

# [fix(c[i], data.c[i]; force=true) for i  in nm_c]
# [fix(b[i], data.b[i]; force=true) for i  in nm_b]
@NLobjective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 

unset_silent(m)
JuMP.optimize!(m)

st0 = termination_status(m)
obj0 = objective_value(m)

psol0 = JuMP.value.(p)
csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
zsol0 = [Asol0...;bsol0...;csol0...;psol0...]
Zsol0 = zsol0*zsol0'
rank(Zsol0)
eigvals(Zsol0)

tol = 1e-5
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol0.<=1)
    @test all(0 .<=bsol0.<=data.p_bar)
    @test all(0 .<=csol0.<=data.assets)
    @test all(0 .<=psol0.<=data.p_bar)
    
    @test norm(min.(p_bar,max.((1+g0).*(Asol0'*psol0[:,1] .+ csol0 .- x0[:,1].*csol0 ) .- g0.*p_bar,0))-psol0[:,1])<tol
end

Id = Matrix{Float64}(I,N,N)
p_lb = min.(p_bar,(Id-Diagonal(x0[:,1]))*(Id-Asol0)\csol0)
L0 = sum(sum( x0[i,j]*csol0[i]+p_bar[i]-psol0[i,j] for i=1:N)/D for j=1:D)
L_lb = sum(sum( x0[i,j]*csol0[i]+p_bar[i]-p_lb[i] for i=1:N)/D for j=1:D)
L_d = sum(x0[:,1]+max.(0,x0[:,1].-w))
@show (L0,L_lb,L_d)

β = (p_bar.-bsol0)./p_bar # financial connectivity: proportion of liabilities inside the network
β⁺ = maximum(β)


expectation(x->x^2, Beta(1.0,2.0))

function clearing_vec!(F, p, x, A)
    F .= p.-min.(p_bar,max.((1+g0)*(A'*p .+csol0 .- x.*csol0) .-g0.*p_bar,0))
end

function shock()
    # draw shocks that produce the probability of default δ 
    a = 1
    b = log.(δ)./(log.(1 .-w./c))
    dist = Beta.(a,[b...])
    draws = rand.(dist, 1)
    vcat(draws...).*c
    #rand(Float64,N) # normal
end

function loss(x,A)
     # total loss in actual network for shock x
    p = nlsolve((F, p)->clearing_vec!(F, p, x, A),[p_bar...], autodiff = :forward)
    @test converged(p)
    sum(x.+p_bar.-p.zero) 
end

function loss_d(x)
    # total loss in disconnected network for shock x
    #b_d = b .+ max.(0, c.-b.-w ) # add fictitious outside liabilities 
    #_d = c .+ max.(0, w.-c.+b ) # add fictitious outside assets
    sum(x+max.(0,x.-w))
end


a = data.assets # total assets
d=  assets .- csol0;# inside assets
f = p_bar .- bsol0;# inside liabilities


# example from Fernando's presentation and Glasserman-Young
using NLsolve
x1 = ones(N)/2
loss(x1,Asol0)
loss_d(x1)

ratio = loss(x1,A0)/loss_d(x1)
bound = 1 + sum(data.delta.*csol0)/((1-β⁺)*sum(csol0))


############### max relaxation ######################
m = Model(with_optimizer(Ipopt.Optimizer,max_iter=10000))

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = data_nm.c[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = data_nm.b[i])   

@constraint(m, sum(A,dims=2).*data.p_bar .== data.p_bar .- b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d

# liabilities are net liabilities: A[i,i]=0 and A[i,j]A[j,i]=0
@constraint(m, [i = 1:N], A[i,i]==0)
for i=1:N
    j=1
    while j < i
        @complements(m, 0 <= A[i,j],  A[j,i] >= 0)
        j += 1
    end
end

# register max and min non-linear functions
maxfun(n1, n2) = max(n1, n2)
minfun(n1, n2) = min(n1, n2)

JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# clearing vector
myexpr = []
for j=1:D
    push!(myexpr, (1+g0).*(A'*p[:,j] .+ c .- x0[:,j].*c ) .- g0.*data.p_bar)
end
    
@variable(m, aux[i=1:N,j=1:D])
for i = 1:N
    for j=1:D
        @constraint(m, aux[i,j] .== myexpr[j][i] )
    end
end

for i = 1:N
    for j=1:D
        @NLconstraint(m,  minfun(data.p_bar[i], aux[i,j]) == p[i,j] ) 
    end
end

# [fix(c[i], data.c[i]; force=true) for i  in nm_c]
# [fix(b[i], data.b[i]; force=true) for i  in nm_b]

#@NLobjective(m, Max , sum(sum( x0[i,j]*c[i]+data.p_bar[i]-p[i,j] for i=1:N)/D for j=1:D) ) #*sum(x[i]+p_bar[i]-p[i] for i=1:N) 
@NLobjective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 

#unset_silent(m)
JuMP.optimize!(m)

termination_status(m)
objective_value(m)



################################
x =  [0.4786082776942293
 0.6261009289132826
 0.0643506085692449
 0.32494449525126634
 0.5115660459579572];
T= Float64
D=1
M = N^2+3*N
γ = 0.0
rng =  Random.seed!(123)
A0 = Asol0

li = LinearIndices((N,N))
liᵀ= transpose(li)
id = li[tril!(trues(N,N), 0)]
idᵀ = liᵀ[tril!(trues(N,N), 0)]
H₃ = sparse(vcat(id,idᵀ),vcat(idᵀ,id),0.5,M,M)
h₃ = spzeros(M)
q₃ = 0.0
h₀ᵀ = sparse(hcat(zeros(1,N^2),zeros(1,N),zeros(1,N),ones(1,N)))
H₄⁽ⁱ⁾ = []
h⁽ⁱ⁾ᵀ = []
q₄⁽ⁱ⁾ = []
h₄⁽ⁱ⁾ᵀ= []
h₄₄⁽ⁱ⁾ᵀ= []
for i=1:N
    push!(H₄⁽ⁱ⁾,sparse( vcat(N*(i-1)+1:N*i,M-N+1:M), vcat(M-N+1:M,N*(i-1)+1:N*i), -(1+γ)*ones(T,2*N)/2,M,M))
    #h₄⁽ⁱ⁾ᵀ = sparse([1],[N^2+N+i],(1+γ)*(1-x[i]),1,M)
    #h₄₄⁽ⁱ⁾ᵀ = sparse([1],[N^2+2*N+i],1.0,1,M)
    h₄⁽ⁱ⁾ᵀ= push!(h₄⁽ⁱ⁾ᵀ,sparse([1],[N^2+N+i],(1+γ)*(1-x[i]),1,M))
    h₄₄⁽ⁱ⁾ᵀ = push!(h₄₄⁽ⁱ⁾ᵀ,sparse([1],[N^2+2*N+i],1.0,1,M))
    push!(h⁽ⁱ⁾ᵀ,-h₄⁽ⁱ⁾ᵀ[i]-h₄₄⁽ⁱ⁾ᵀ[i])
    push!(q₄⁽ⁱ⁾,γ*p_bar[i])
end

tol = 1e-7
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol0.<=1)
    @test all(0 .<=bsol0.<=data.p_bar)
    @test all(0 .<=csol0.<=data.assets)
    @test all(0 .<=psol0.<=data.p_bar)
    @test norm(min.(data.p_bar,max.((1+g0).*(Asol0'*psol0[:,1] .+ csol0 .- x0[:,1].*csol0 ) .- g0.*data.p_bar,0))-psol0[:,1])<tol
    
    @test all(eigvals(Zsol0).>-tol)
    @test all(eigvals([Zsol0 zsol0;zsol0' 1.0]).>-tol)
    @test norm(dot(H₃, Zsol0))<tol
    
    @test (dot(H₄⁽ⁱ⁾[1], Zsol0)+(h⁽ⁱ⁾ᵀ[1]*zsol0)[1] - q₄⁽ⁱ⁾[1])<=0
    @test (dot(H₄⁽ⁱ⁾[2], Zsol0)+(h⁽ⁱ⁾ᵀ[2]*zsol0)[1] - q₄⁽ⁱ⁾[2])<=0
    @test (dot(H₄⁽ⁱ⁾[3], Zsol0)+(h⁽ⁱ⁾ᵀ[3]*zsol0)[1] - q₄⁽ⁱ⁾[3])<=0
    @test (dot(H₄⁽ⁱ⁾[4], Zsol0)+(h⁽ⁱ⁾ᵀ[4]*zsol0)[1] - q₄⁽ⁱ⁾[4])<=0
    @test (dot(H₄⁽ⁱ⁾[5], Zsol0)+(h⁽ⁱ⁾ᵀ[5]*zsol0)[1] - q₄⁽ⁱ⁾[5])<=0
end

m = Model()
set_optimizer(m, CSDP.Optimizer) # SCS, ProxSDP, CSDP
set_optimizer_attribute(m, MOI.Silent(),false)

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = psol0[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = csol0[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = bsol0[i])  
@variable(m, Z[i=1:N^2+3*N,j=1:N^2+3*N], PSD, start = Zsol0[i,j])
z = vcat(A...,b,c,p)

@constraint(m, sum(A,dims=2).*data.p_bar .== data.p_bar .- b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d
@constraint(m, dot(H₃, Z) == 0.0 )
# for i = 1:N
#     @constraint(m, dot(H₄⁽ⁱ⁾[i], Z)+(h⁽ⁱ⁾ᵀ[i]*z)[1] <= q₄⁽ⁱ⁾[i])
# end
@constraint(m, dot(H₄⁽ⁱ⁾[1], Z)+(h⁽ⁱ⁾ᵀ[1]*z)[1] - q₄⁽ⁱ⁾[1] <= 0.0)
@constraint(m, dot(H₄⁽ⁱ⁾[2], Z)+(h⁽ⁱ⁾ᵀ[2]*z)[1] - q₄⁽ⁱ⁾[2] <= 0.0)
@constraint(m, dot(H₄⁽ⁱ⁾[3], Z)+(h⁽ⁱ⁾ᵀ[3]*z)[1] - q₄⁽ⁱ⁾[3] <= 0.0)
@constraint(m, dot(H₄⁽ⁱ⁾[4], Z)+(h⁽ⁱ⁾ᵀ[4]*z)[1] - q₄⁽ⁱ⁾[4] <= 0.0)
@constraint(m, dot(H₄⁽ⁱ⁾[5], Z)+(h⁽ⁱ⁾ᵀ[5]*z)[1] - q₄⁽ⁱ⁾[5] <= 0.0)

@constraint(m, Symmetric([Z z;z' 1.0]) in PSDCone())

#@objective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 
@objective(m, Min , sum(z[end-N+1:end]) )
#@objective(m, MathOptInterface.FEASIBILITY_SENSE ) 

JuMP.optimize!(m)

termination_status(m)
objective_value(m)

psolJ = JuMP.value.(p)
csolJ = JuMP.value.(c)
bsolJ = JuMP.value.(b)
AsolJ = JuMP.value.(A)
zsolJ = JuMP.value.(z)
ZsolJ = JuMP.value.(Z)
zsolJ*zsolJ' .- ZsolJ
rank(ZsolJ)
eigvals(ZsolJ)


tol = 1e-4
@testset "check solution" begin
    @test norm( sum(AsolJ,dims=2).* data.p_bar .- (data.p_bar .- bsolJ)) < tol
    @test norm( AsolJ' * data.p_bar .- (data.assets .- csolJ)) < tol
    @test norm(diag(AsolJ)) < tol
    @test norm([AsolJ[i,j]*AsolJ[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=AsolJ.<=1)
    @test all(0 .<=bsolJ.<=data.p_bar)
    @test all(0 .<=csolJ.<=data.assets)
    @test all(0 .<=psolJ.<=data.p_bar)
    @test norm(min.(data.p_bar,max.((1+γ).*(AsolJ'*psolJ[:,1] .+ csolJ .- x0[:,1].*csolJ ) .- γ.*data.p_bar,0))-psolJ[:,1])<tol
    
    @test all(eigvals(ZsolJ).>-tol)
    @test all(eigvals([ZsolJ zsolJ;zsolJ' 1.0]).>-tol)
    @test norm(dot(H₃, ZsolJ))<tol
    
    @test (dot(H₄⁽ⁱ⁾[1], ZsolJ)+(h⁽ⁱ⁾ᵀ[1]*zsolJ)[1] - q₄⁽ⁱ⁾[1])<=0
    @test (dot(H₄⁽ⁱ⁾[2], ZsolJ)+(h⁽ⁱ⁾ᵀ[2]*zsolJ)[1] - q₄⁽ⁱ⁾[2])<=0
    @test (dot(H₄⁽ⁱ⁾[3], ZsolJ)+(h⁽ⁱ⁾ᵀ[3]*zsolJ)[1] - q₄⁽ⁱ⁾[3])<=0
    @test (dot(H₄⁽ⁱ⁾[4], ZsolJ)+(h⁽ⁱ⁾ᵀ[4]*zsolJ)[1] - q₄⁽ⁱ⁾[4])<=0
    @test (dot(H₄⁽ⁱ⁾[5], ZsolJ)+(h⁽ⁱ⁾ᵀ[5]*zsolJ)[1] - q₄⁽ⁱ⁾[5])<=0
end


################################

γ = 0
T = Float64
N = length(p_bar)
x =  [0.4786082776942293
 0.6261009289132826
 0.0643506085692449
 0.32494449525126634
 0.5115660459579572];
M = N^2+3*N

low = zeros(T,M)
upp = vcat(ones(T,N^2),p_bar,assets,p_bar)

H₁ = spzeros(M,M)
h₁₁ᵀ = hcat(kron(Matrix{T}(I,N,N),p_bar'))
h₁ᵀ = hcat(h₁₁ᵀ,zeros(T,N,N),Matrix{T}(I,N,N),zeros(T,N,N))
q₁ = assets

H₂ = spzeros(M,M)
h₂₁ᵀ = repeat(diagm(p_bar),1,N)
h₂ᵀ = hcat(h₂₁ᵀ,Matrix{T}(I,N,N),zeros(T,N,N),zeros(T,N,N))
q₂ = p_bar

li = LinearIndices((N,N))
liᵀ= transpose(li)
id = li[tril!(trues(N,N), 0)]
idᵀ = liᵀ[tril!(trues(N,N), 0)]
H₃= sparse(vcat(id,idᵀ),vcat(idᵀ,id),0.5,M,M)
h₃ = zeros(M)
q₃ = 0.0

H₄⁽ⁱ⁾ = []
h⁽ⁱ⁾ᵀ = []
q₄⁽ⁱ⁾ = []
for i=1:N
    push!(H₄⁽ⁱ⁾,sparse( vcat(N*(i-1)+1:N*i,M-N+1:M), vcat(M-N+1:M,N*(i-1)+1:N*i), -(1+γ)*ones(T,2*N)/2,M,M))
    h₄⁽ⁱ⁾ᵀ = sparse([1],[N^2+N+i],(1+γ)*(1-x[i]),1,M)
    h₄₄⁽ⁱ⁾ᵀ = sparse([1],[N^2+2*N+i],1.0,1,M)
    push!(h⁽ⁱ⁾ᵀ,-h₄⁽ⁱ⁾ᵀ-h₄₄⁽ⁱ⁾ᵀ)
    push!(q₄⁽ⁱ⁾,γ*p_bar[i])
end

h₀ᵀ = hcat(zeros(T,1,N^2),zeros(T,1,N),zeros(T,1,N),ones(T,1,N))

# spy(H₀)
# spy(H₀[1:N^2,1:N^2])

z = Variable(M,Positive())
Z = Variable(M,M,Positive())

problem = minimize(h₀ᵀ*z)
problem.constraints = [
                        z>=low,z<=upp,
                        h₁ᵀ*z == q₁,
                        h₂ᵀ*z == q₂,
                        tr(H₃*Z) == q₃,
                        [tr(H₄⁽ⁱ⁾[i]*Z)+(h⁽ⁱ⁾ᵀ[i]*z)[1] <= q₄⁽ⁱ⁾[i] for i in 1:N]...,
                        [Z z; z' 1.0] ⪰ 0.0,
                        abs(eigmax([Z z; z' 1.0]))<=1.1
                       ]
set_value!(Z,zsol0*zsol0')
set_value!(z,zsol0)
Convex.solve!(problem,  SCS.Optimizer(verbose=true, eps=1e-10,linear_solver=SCS.DirectSolver), warmstart=true)
Convex.solve!(problem,  ProxSDP.Optimizer())
Convex.solve!(problem,  CSDP.Optimizer())
Convex.solve!(problem,  Tulip.Optimizer())

zsol = evaluate(z)
Zsol = evaluate(Z)

Asol = reshape(zsol[1:N^2],N,N)
bsol = zsol[N^2+1:N^2+N]
csol = zsol[N^2+N+1:N^2+2*N]
psol = zsol[N^2+2*N+1:N^2+3*N]
max(abs.(zsol*zsol' .- Zsol)...)
rank([Zsol zsol; zsol' 1.0])

eigvals([Zsol zsol; zsol' 1.0])


tol = 1e-4
@testset "check solution" begin
    @test norm( sum(Asol,dims=2).* data.p_bar .- (data.p_bar .- bsol)) < tol
    @test norm( Asol' * data.p_bar .- (data.assets .- csol)) < tol
    @test norm(diag(Asol)) < tol
    @test norm([Asol[i,j]*AsolJ[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol.<=1)
    @test all(0 .<=bsol.<=data.p_bar)
    @test all(0 .<=csol.<=data.assets)
    @test all(0 .<=psol.<=data.p_bar)
    @test norm(min.(data.p_bar,max.((1+γ).*(Asol'*psol[:,1] .+ csol .- x0[:,1].*csol ) .- γ.*data.p_bar,0))-psol[:,1])<tol
    
    @test all(real(eigvals(Zsol)).>-tol)
    @test norm(imag(eigvals(Zsol)))<tol
    @test all(eigvals([Zsol zsol;zsol' 1.0]).>-tol)
    @test norm(dot(H₃, Zsol))<tol
    
    @test (dot(H₄⁽ⁱ⁾[1], Zsol)+(h⁽ⁱ⁾ᵀ[1]*zsol)[1] - q₄⁽ⁱ⁾[1])<=0
    @test (dot(H₄⁽ⁱ⁾[2], Zsol)+(h⁽ⁱ⁾ᵀ[2]*zsol)[1] - q₄⁽ⁱ⁾[2])<=0
    @test (dot(H₄⁽ⁱ⁾[3], Zsol)+(h⁽ⁱ⁾ᵀ[3]*zsol)[1] - q₄⁽ⁱ⁾[3])<=0
    @test (dot(H₄⁽ⁱ⁾[4], Zsol)+(h⁽ⁱ⁾ᵀ[4]*zsol)[1] - q₄⁽ⁱ⁾[4])<=0
    @test (dot(H₄⁽ⁱ⁾[5], Zsol)+(h⁽ⁱ⁾ᵀ[5]*zsol)[1] - q₄⁽ⁱ⁾[5])<=0
end


##
Z = Variable(M+1,M+1,Positive())
r = Variable(1,Positive())
lowZ = [low.^2; 1.0]
uppZ = [upp.^2; 1.0]
hh₀ᵀ = [zeros(M,M) h₀ᵀ'/2; h₀ᵀ/2 0.0]
HH₃ = [H₃ zeros(M); zeros(1,M) 0.0]
hh₁ᵀ = [];
hh₂ᵀ = [];
HH₄⁽ⁱ⁾ = [];
for i =1:N
    push!(hh₁ᵀ, [zeros(M,M) view(h₁ᵀ,i,:)/2; view(h₁ᵀ,i,:)'/2 0.0] )
    push!(hh₂ᵀ, [zeros(M,M) view(h₂ᵀ,i,:)/2; view(h₂ᵀ,i,:)'/2 0.0] )
    push!(HH₄⁽ⁱ⁾ , [H₄⁽ⁱ⁾[i] (1/2)*h⁽ⁱ⁾ᵀ[i]'; (1/2)*h⁽ⁱ⁾ᵀ[i] 0.0])
end

#problem = minimize(tr(hh₀ᵀ*Z))
#problem = minimize(r)
problem = minimize(r)
problem.constraints = [
                        diag(Z)>=lowZ,diag(Z)<=uppZ,
                        [tr(hh₁ᵀ[i]*Z) <= q₁[i] for i in 1:N]...,
                        [tr(hh₂ᵀ[i]*Z) <= q₂[i] for i in 1:N]...,
                        tr(HH₃*Z) == q₃,
                        [tr(HH₄⁽ⁱ⁾[i]*Z) <= q₄⁽ⁱ⁾[i] for i in 1:N]...,
                        Z ⪰ 0.0,
                        #nuclearnorm(Z)<=r
                       ]
set_value!(r,4);
set_value!(Z,vcat(zsol0,1.0)*vcat(zsol0,1.0)');
Convex.solve!(problem,  SCS.Optimizer(verbose=true, eps=1e-10,linear_solver=SCS.DirectSolver), warmstart=true, verbose=true)
Convex.solve!(problem,  COSMO.Optimizer(eps_rel=1e-10,eps_abs=1e-10), warmstart=true, verbose=true)
# Convex.solve!(problem,  COSMO.Optimizer(eps_rel=1e-10,eps_abs=1e-10,merge_strategy = COSMO.CliqueGraphMerge,decompose = true, kkt_solver = with_options(CholmodKKTSolver)), warmstart=true, verbose=true)
Convex.solve!(problem,  ProxSDP.Optimizer(),  verbose=true)


rsolH = evaluate(r)
ZsolH = evaluate(Z)
rank(ZsolH[1:M,1:M])
eigvals(ZsolH[1:M,1:M])
@show  "Rank", sum(abs.(eigvals(ZsolH[1:M,1:M])).>1e-5)
eigZH = eigen(ZsolH)
eigZH.values[end]
#ZeigH = eigZH.vectors[:,end]
ZeigH = sqrt.(max.(0.0,diag(ZsolH)))
AsolH = real(reshape(ZeigH[1:N^2],N,N))
bsolH = real(ZeigH[N^2+1:N^2+N])
csolH = real(ZeigH[N^2+N+1:N^2+2*N])
psolH = real(ZeigH[N^2+2*N+1:N^2+3*N])
tr(hh₀ᵀ*ZsolH)
sum(psolH)

tol = 1e-4
@testset "check solution" begin
    @test norm( sum(Asol,dims=2).* data.p_bar .- (data.p_bar .- bsol)) < tol
    @test norm( Asol' * data.p_bar .- (data.assets .- csol)) < tol
    @test norm(diag(Asol)) < tol
    @test norm([Asol[i,j]*AsolJ[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol.<=1)
    @test all(0 .<=bsol.<=data.p_bar)
    @test all(0 .<=csol.<=data.assets)
    @test all(0 .<=psol.<=data.p_bar)
    @test norm(min.(data.p_bar,max.((1+γ).*(Asol'*psol[:,1] .+ csol .- x0[:,1].*csol ) .- γ.*data.p_bar,0))-psol[:,1])<tol
    
    @test all(real(eigvals(Zsol)).>-tol)
    @test norm(imag(eigvals(Zsol)))<tol
    @test all(eigvals([Zsol zsol;zsol' 1.0]).>-tol)
    @test norm(dot(H₃, Zsol))<tol
    
    @test (dot(H₄⁽ⁱ⁾[1], Zsol)+(h⁽ⁱ⁾ᵀ[1]*zsol)[1] - q₄⁽ⁱ⁾[1])<=0
    @test (dot(H₄⁽ⁱ⁾[2], Zsol)+(h⁽ⁱ⁾ᵀ[2]*zsol)[1] - q₄⁽ⁱ⁾[2])<=0
    @test (dot(H₄⁽ⁱ⁾[3], Zsol)+(h⁽ⁱ⁾ᵀ[3]*zsol)[1] - q₄⁽ⁱ⁾[3])<=0
    @test (dot(H₄⁽ⁱ⁾[4], Zsol)+(h⁽ⁱ⁾ᵀ[4]*zsol)[1] - q₄⁽ⁱ⁾[4])<=0
    @test (dot(H₄⁽ⁱ⁾[5], Zsol)+(h⁽ⁱ⁾ᵀ[5]*zsol)[1] - q₄⁽ⁱ⁾[5])<=0
    
    sum(eigvals(ZsolH))
end

@testset "does true solution satisfy relaxation?" begin
    zsol00 = vcat(zsol0,1.0)
    Zsol00 = zsol00*zsol00';
    F = svd(Zsol00)
    @test all(diag(Zsol00).>=lowZ)
    @test all(diag(Zsol00).<=uppZ)
    @test norm([tr(hh₁ᵀ[i]*Zsol00)-q₁[i] for i in 1:N])<tol
    @test norm([tr(hh₂ᵀ[i]*Zsol00)-q₂[i] for i in 1:N])<tol
    @test abs(tr(HH₃*Zsol00)-q₃)<tol
    @test all([tr(HH₄⁽ⁱ⁾[i]*Zsol00)-q₄⁽ⁱ⁾[i] for i in 1:N].<0)
    @test all(real(eigvals(Zsol00)).>-tol)
    @show sum(eigvals(Zsol00))
    @show sum(F.S)
end
#########################################
#fix!(x, v),free!(x)


###############################################
using ProximalAlgorithms
using ProximalOperators

  using RecursiveArrayTools: ArrayPartition

diag(Z)>=lowZ,diag(Z)<=uppZ,
[tr(hh₁ᵀ[i]*Z) <= q₁[i] for i in 1:N]...,
[tr(hh₂ᵀ[i]*Z) <= q₂[i] for i in 1:N]...,
tr(HH₃*Z) == q₃,
[tr(HH₄⁽ⁱ⁾[i]*Z) <= q₄⁽ⁱ⁾[i] for i in 1:N]...,
Z ⪰ 0.0

M = N^2 + 3*N
tol=1e-8
maxit=50000
  A = rand(N,N);
  b = rand(N);
  c = rand(N);
  p = rand(N);
  x = (A,b,c,p)
  u = zero(x)
  z = deepcopy(x)
  f = NormL1(1.0)
  g = IndBallRank(1)
  gam = 0.1
  for it = 1:maxit
    # perform f-update step
    prox!(p, f, z[4] - u[4], gam)
    
    # perform g-update step
    prox!(A, g, x[1] + u[1], gam)
     x = (A,b,c,p);
    # stopping criterion
    if norm(x-z, Inf) <= tol*(1+norm(u, Inf))
      break
    end
    # dual update
    u .+= x - z
  end
  return z
end
m, n, k, sig = 500, 2500, 100, 1e-3
A = randn(m, n)
x_true = [randn(k)..., zeros(n-k)...]
b = A*x_true + sig*randn(m)
lam = 0.1*norm(A'*b, Inf)
x_admm = lasso_admm(A, b, lam, zeros(n))


f = SeparableSum(NormL1(), NuclearNorm());
x = randn(10); # some random vector
Y = randn(20, 30); # some random matrix
f_xY = f((x, Y)) # evaluates f at (x, Y)
(u, V), f_uV = prox(f, (x, Y), 1.3); # computes prox at (x, Y)

###################################


constraint8 = BilinearConstraint(vec(A)',H,vec(A),0.0,X=-zeros(1,N^2)/2,Y=-zeros(N^2,1)/2,λ=15.0, W1=Matrix(I,1,1)/2, W2=Matrix(I,1,1)/2)
bp=BilinearProblem(problem,constraint8)
SequentialConvexRelaxation.solve!(bp,()->SCS.Optimizer(verbose=1,acceleration_lookback=1,max_iters=100000,eps=1e-8),iterations=10,update_weights=true,weight_update_tuning_param=0.1)
#constraint9 = BinaryConstraint(aux,X=0.4,λ = 2.0)

SequentialConvexRelaxation.solve!(bp,()->SCS.Optimizer,iterations=10,update_weights=true,weight_update_tuning_param=0.1)


Asol = evaluate(A)
vec(Asol)'*H*vec(Asol)
tol = 1e-6
@testset "convex_problem" begin
    @test all(Asol.-1 .<= tol)
    @test all(Asol .>= -tol)
    @test theta*auxsol-Asol[1,2]>=-tol
    @test theta*(1-auxsol)-Asol[2,1]>=-tol
end


Aplot = copy(Asol)
Aplot[Aplot.<1f-4].=0
Aplot


julia> Pkg.clone("https://github.com/rdeits/CouenneNL.jl")
julia> Pkg.build("CouenneNL")
PiecewiseLinearOpt.jl




import Pkg
Pkg.add("Pajarito")
Pkg.add("GLPKMathProgInterface")
Pkg.add("Convex")
using Pajarito, GLPKMathProgInterface, SCS #load packages
using GLPK
mysolver = PajaritoSolver(log_level = 3, #use verbose output
mip_solver = GLPKSolverMIP(msg_lev = GLPK.MSG_OFF), #set MIP solver
cont_solver = SCSSolver(eps = 1e-6, verbose = 0)) #set conic solver


using Convex, LinearAlgebra
p=1;m=2;n=3;V=rand(n,p)
mp = Variable(p, Positive(), :Int) #create p nonneg. integer variables
eOpt = maximize(eigmin(V * diag(mp./m) * V'), #max. min. eigenvalue
sum(mp) <= m) #add linear constraint
solve!(eOpt, Pajarito.mysolver()) #solve model using Pajarito solver
@show eOpt.status, eOpt.optval, mp.value #show solve status and results



Pkg.add("Pavito")
using Convex, Pavito

# Set up Convex.jl model, solve, print solution
function gatesizing(yUB, solver)
    fe = [1, 0.8, 1, 0.7, 0.7, 0.5, 0.5] .* [1, 2, 1, 1.5, 1.5, 1, 2]
    Cout6 = 10
    Cout7 = 10

    y = Variable(7, Positive())
    z = Variable(yUB, 7, Positive(), :Bin)

    D1 = exp(-y[1]) + exp(-y[1] + y[4])
    D2 = 2 * exp(-y[2]) + exp(-y[2] + y[4]) + exp(-y[2] + y[5])
    D3 = 2 * exp(-y[3]) + exp(-y[3] + y[5]) + exp(-y[3] + y[7])
    D4 = 2 * exp(-y[4]) + exp(-y[4] + y[6]) + exp(-y[4] + y[7])
    D5 = exp(-y[5]) + exp(-y[5] + y[7])
    D6 = Cout6 * exp(-y[6])
    D7 = Cout7 * exp(-y[7])

    P = minimize(
        maximum([(D1+D4+D6), (D1+D4+D7), (D2+D4+D6), (D2+D4+D7), (D2+D5+D7), (D3+D5+D6), (D3+D7)]),
        sum(fe .* exp(y)) <= 20,
        sum(exp(y)) <= 100)

    for i in 1:7
        P.constraints += (sum(z[:,i]) == 1)
        P.constraints += (y[i] == sum([log(j) * z[j,i] for j=1:yUB]))
    end

    solve!(P, solver)

    println("\nCircuit delay (obj) = $(P.optval)")
    println("Scale factors (exp(y)):\n$(exp.(y.value))")
    println("Value indicators (z):\n$(round.(z.value))")
end


#=========================================================
Choose solvers and options
=========================================================#

mip_solver_drives = false
rel_gap = 1e-5


using Cbc
mip_solver = CbcSolver()

using CPLEX
mip_solver = CplexSolver(
    CPX_PARAM_SCRIND=(mip_solver_drives ? 1 : 0),
    CPX_PARAM_EPINT=1e-8,
    CPX_PARAM_EPRHS=1e-7,
    CPX_PARAM_EPGAP=(mip_solver_drives ? 1e-5 : 1e-9)
)


using Ipopt
cont_solver = IpoptSolver(print_level=0)

solver = PavitoSolver(
    mip_solver_drives=mip_solver_drives,
    log_level=2,
    rel_gap=rel_gap,
	mip_solver=mip_solver,
	cont_solver=cont_solver,
)


#=========================================================
Specify data
=========================================================#

# Integer upper bound on integer gate size factors
# Optimal solution for this instances has all exp(y) <= 3,
# so yUB larger than 2 is not constraining, but impacts number of variables
yUB = 3


#=========================================================
Solve Convex.jl model
=========================================================#

# Optimal solution is [2,3,3,3,2,3,3] with value 8.3333
gatesizing(yUB, solver)





function lasso_admm(A, b, lam, x; tol=1e-8, maxit=50000)
  u = zero(x)
  z = copy(x)
  f = LeastSquares(A, b)
  g = NormL1(lam)
  gam = 100.0/norm(A)^2
  for it = 1:maxit
    # perform f-update step
    prox!(x, f, z - u, gam)
    # perform g-update step
    prox!(z, g, x + u, gam)
    # stopping criterion
    if norm(x-z, Inf) <= tol*(1+norm(u, Inf))
      break
    end
    # dual update
    u .+= x - z
  end
  return z
end


import Pkg
Pkg.add("ProximalAlgorithms")
Pkg.add("ProximalOperators")
using ProximalAlgorithms
using ProximalOperators
using LinearAlgebra
using Random
using Test

@testset "Nonconvex QP (tiny, $T)" for T in [Float64]
    Q = Matrix(Diagonal(T[-0.5, 1.0]))
    q = T[0.3, 0.5]
    low = T(-1.0)
    upp = T(+1.0)

    f = Quadratic(Q, q)
    g = IndBox(low, upp)

    n = 2

    Lip = maximum(diag(Q))
    gamma = T(0.95) / Lip

    @testset "PANOC" begin
        x0 = zeros(T, n)
        solver = ProximalAlgorithms.PANOC{T}()
        x, it = solver(x0, f = f, g = g)
        z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
        @test norm(x - z, Inf) / gamma <= solver.tol
    end

    @testset "ZeroFPR" begin
        x0 = zeros(T, n)
        solver = ProximalAlgorithms.ZeroFPR{T}()
        x, it = solver(x0, f = f, g = g)
        z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
        @test norm(x - z, Inf) / gamma <= solver.tol
    end

    @testset "LiLin" begin
        x0 = zeros(T, n)
        solver = ProximalAlgorithms.LiLin{T}(gamma = gamma)
        x, it = solver(x0, f = f, g = g)
        z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
        @test norm(x - z, Inf) / gamma <= solver.tol
    end
end

@testset "Nonconvex QP (small, $T)" for T in [Float64]
    @testset "Random problem $k" for k = 1:5
        Random.seed!(k)

        n = 100
        A = randn(T, n, n)
        U, R = qr(A)
        eigenvalues = T(2) .* rand(T, n) .- T(1)
        Q = U * Diagonal(eigenvalues) * U'
        Q = 0.5 * (Q + Q')
        q = randn(T, n)

        low = T(-1.0)
        upp = T(+1.0)

        f = Quadratic(Q, q)
        g = IndBox(low, upp)

        Lip = maximum(abs.(eigenvalues))
        gamma = T(0.95) / Lip

        TOL = 1e-4

        @testset "PANOC" begin
            x0 = zeros(T, n)
            solver = ProximalAlgorithms.PANOC{T}(tol = TOL)
            x, it = solver(x0, f = f, g = g)
            z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
            @test norm(x - z, Inf) / gamma <= solver.tol
        end

        @testset "ZeroFPR" begin
            x0 = zeros(T, n)
            solver = ProximalAlgorithms.ZeroFPR{T}(tol = TOL)
            x, it = solver(x0, f = f, g = g)
            z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
            @test norm(x - z, Inf) / gamma <= solver.tol
        end

        @testset "LiLin" begin
            x0 = zeros(T, n)
            solver = ProximalAlgorithms.LiLin(gamma = gamma, tol = TOL)
            x, it = solver(x0, f = f, g = g)
            z = min.(upp, max.(low, x .- gamma .* (Q * x + q)))
            @test norm(x - z, Inf) / gamma <= solver.tol
        end
    end
end


indA = 1:N^2
indb = N^2+1:N^2+N
indc = N^2+N+1:N^2+2*N
indp = N^2+2*N+1:N^2+3*N
ind_ASOC = N^2+3*N+1:N^2+3*N+2
M = ind_ASOC[end]
id = Matrix{Float64}(I,M,M)

lb = COSMO.Constraint(id, zeros(M), COSMO.Nonnegatives)
ub = COSMO.Constraint(Matrix{Float64}(I,M-2,M-2), -vcat(ones(indA[end]),p_bar,assets,p_bar), COSMO.Nonnegatives,M,1:M-2)

aux1 = hcat(repeat(diagm(p_bar),outer=[1,N]),Matrix{Float64}(I,N,N))
aux2 = hcat(repeat(diagm(p_bar),inner=[1,N]),zeros(N,N),Matrix{Float64}(I,N,N))

cons_liab = COSMO.Constraint(aux1, -p_bar, COSMO.ZeroSet,M,indA[1]:indb[end])
cons_assets = COSMO.Constraint(aux2, -assets, COSMO.ZeroSet,M,indA[1]:indc[end])

aux3 = zeros(2,M)
aux3[1,end-1]=1
aux3[2,end]=1
cons_q = COSMO.Constraint(aux3, zeros(2), COSMO.SecondOrderCone(2))

aux4 = zeros(2,M)
aux4[1,4]=1;aux4[1,end-1]=-1/sqrt(2);aux4[1,end]=-1/sqrt(2)
aux4[1,2]=1;aux4[1,end-1]=-1/sqrt(2);aux4[1,end]=1/sqrt(2)
cons_q2 = COSMO.Constraint(aux4, zeros(2), COSMO.ZeroSet)

constraints = [lb; ub; cons_liab; cons_assets; cons_q; cons_q2]

Linear(A,b)


###################################
import Pkg; Pkg.add("SumOfSquares"); Pkg.add("DynamicPolynomials"); Pkg.add("PolyJuMP")
using PolyJuMP
using SumOfSquares
using DynamicPolynomials
using ProxSDP, SCS

# Using ProxSDP as the SDP solver
model = SOSModel(with_optimizer(ProxSDP.Optimizer))
model = SOSModel(with_optimizer(SCS.Optimizer, max_iters=100000))
set_optimizer_attribute(model, MOI.Silent(),false)

@polyvar x z
@variable(model, t)
p = x^4 + x^2 - 3*x^2*z^2 + z^6 - t
@constraint(model, x*z == 0)
@constraint(model, p >= 0)
@objective(model, Max, t)
optimize!(model)
println("Solution: $(value(t))")
value(x)


m = SOSModel(with_optimizer(ProxSDP.Optimizer))
set_optimizer_attribute(m, MOI.Silent(),false)
@polyvar a12 a21 a31 a13
X = monomials([x, y], 0:2)
@variable(model, p, Poly(X))
poly = a12*a21+a13*a31
@constraint(m,poly==0.0)
@constraint(m,a12<=1.0)
@constraint(m,a21<=1.0)
@constraint(m,a13<=1.0)
@constraint(m,a31<=1.0)
@constraint(m,a12>=0.0)
@constraint(m,a21>=0.0)
@constraint(m,a13>=0.0)
@constraint(m,a31>=0.0)
@variable(m, aux)
@constraint(m, aux==a12)
@objective(m, Max , aux) 
unset_silent(m)
JuMP.optimize!(m)











@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = psol0[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=Asol0[i,j])  
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = csol0[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = bsol0[i])   

@constraint(m, sum(A,dims=2).*data.p_bar .== data.p_bar .- b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d

# liabilities are net liabilities: A[i,i]=0 and A[i,j]A[j,i]=0
@constraint(m, [i = 1:N], A[i,i]==0)

# @constraint(m,a12==A[1,2])
# @constraint(m,a21==A[2,1])

# register max and min non-linear functions
maxfun(n1, n2) = max(n1, n2)
minfun(n1, n2) = min(n1, n2)

JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# clearing vector
myexpr = []
for j=1:D
    push!(myexpr, (1+g0).*(A'*p[:,j] .+ c .- x0[:,j].*c ) .- g0.*data.p_bar)
end
    
@variable(m, aux[i=1:N,j=1:D])
for i = 1:N
    for j=1:D
        @constraint(m, aux[i,j] .== myexpr[j][i] )
    end
end

for i = 1:N
    for j=1:D
        @NLconstraint(m,  minfun(data.p_bar[i], maxfun(aux[i,j],0)) == p[i,j] ) 
    end
end

@NLobjective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 



termination_status(m)
objective_value(m)

JuMP.value.(p)
JuMP.value.(c)
JuMP.value.(b)
JuMP.value.(A)
[Asol0...;bsol0...;csol0...;psol0...]
zsol0*zsol0'
rank(Zsol0)
eigvals(Zsol0)




using DynamicPolynomials, SemialgebraicSets
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SemialgebraicSets
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)

model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)

@show termination_status(model)
@show value(a)
@show value(b)
@show objective_value(model)

using CSDP
optimizer = optimizer_with_attributes(CSDP.Optimizer);
model = SOSModel(optimizer)
@polyvar a12 a21 a13 a31
X = monomials([a12,a21,a13,a31], 0:2)
p = a12*a21+a13*a31
S = @set a12 >= 0.0 && a21 >= 0.0 && a13 >= 0.0 && a31 <= 0.0 && a12 <= 1.0 && a21<=1.0 && a31<=1.0  && a13<=1.0
@variable(model, α)
#@variable(model, β)
@objective(model, Max, α )
#@constraint(model, c3,  β==a12, domain = S, maxdegree = 3)
@constraint(model, c4, p >= α, domain = S, maxdegree = 20)
optimize!(model)


@show termination_status(model)
@show objective_value(model)
ν3 = moment_matrix(c3)
extractatoms(ν3, 1e-3) 

