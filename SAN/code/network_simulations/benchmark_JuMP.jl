using LinearAlgebra, DataFrames, XLSX, Missings, JuMP, Ipopt, Random, Complementarity, LinearAlgebra, Test, Expectations, Distributions, NLsolve, DistributionsAD, SpecialFunctions, Quadrature
using Quadrature, ForwardDiff, FiniteDiff, Zygote, Cuba, DistributionsAD, Cubature,HCubature, Test, LinearAlgebra

N = 2 # keep largest `N' nodes by assets

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... )
unique!(data) # delete duplicate rows, use `nonunique(data)' to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
#N = size(data,1) # number of nodes
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
show(data, true)

p_bar = data_nm.p_bar
assets = data_nm.assets
delta = data_nm.delta
w= data_nm.w
g0 = 0.0   # bankruptcy cost

rng =  Random.seed!(123)
# α=1.0f0 ;β=2.0f0 ;
# D = 1 # number of draws
# x0 = rand(N,D)
A0 = zeros(N,N)
c0 =  data_nm.c

# α0= []
# β0= []
# for i=1:N
#      function ff!(F, x)
#         F[1] = 1-cdf(Beta(max(0.1,x[1]),max(0.1,x[2])),w[i]/c0[i])-delta[i] 
#     end
#     sol = nlsolve(ff!,[2.0,30.0]; method = :newton)
#     push!(α0,sol.zero[1])
#     push!(β0,sol.zero[2])
# end
α0 = 1.5*ones(N)
β0= []
for i=1:N
    function ff!(F, x)
        F[1] = 1-cdf(Beta(α0[i],x[1]),w[i]/c0[i])-delta[i] 
    end
    if i==1
        sol = nlsolve(ff!,[30.0])
    else
        sol = nlsolve(ff!,[β0[i-1]])
    end
    push!(β0,sol.zero[1])
end

dist = Beta.(α0,β0);

dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))
function dist_cdf(p...)
    α,β,wc = p[1], p[2], p[3]
    s(x,wc) = 1 ./(1 .+exp.(-1000*(x-wc))) # approximates indicator function of x>=wc
    obj(x,p) = dist_pdf(x,p)*s(x,wc...)
    prob = QuadratureProblem(obj,0.0,1.0,[α...,β...])
    Quadrature.solve(prob,QuadGKJL(),reltol=1e-5,abstol=1e-5)[1]
end

function dist_pdf_mv(x,p...)
    #α,β,x = p[1:N], p[N+1:2*N], p[2*N+1:3*N]
    α,β = p[1][1:N], p[1][N+1:2*N]
    joint_pdf = 1.0;
    for i=1:N
       joint_pdf = joint_pdf*dist_pdf(x[i],[α[i],β[i]]) 
    end
    sum(x)*joint_pdf
end
function ev(p...)
    prob2 = QuadratureProblem(dist_pdf_mv,zeros(N),ones(N),p...)
    Quadrature.solve(prob2,CubaCuhre(),reltol=1e-5,abstol=1e-5)[1]
end

ev([α0...,β0...])
sum(α0./(α0.+β0))

# set up optimization
# m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps",max_iter=10000,hessian_approximation="limited-memory"))
m = Model(with_optimizer(Ipopt.Optimizer,max_iter=100000,print_level=6))
@variable(m, 0.8<=α[1:N]<=3.0, start = 1.5) 
@variable(m, 1.0<=β[1:N]<=150.0, start = 30.0) 
@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, w[i]<=c[i=1:N]<=data.assets[i], start = c0[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = data.p_bar[i]/2)   

@constraint(m, -sum(A,dims=2).*data.p_bar .+ data.p_bar .==  b ) # payments to other nodes add up to inside liabilities f
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

# match probabilities of default
JuMP.register(m, :dist_cdf, 3, dist_cdf, autodiff=true)
@NLexpression(m,wc[i=1:N],w[i]/c[i])
for i in 1:N
    @eval ($(Symbol("vv$i"))) = [α[$i],β[$i],wc[$i]]
    @eval @NLconstraint(m,dist_cdf(($(Symbol("vv$i")))...)== delta[$i]) 
end

# [fix(c[i], data.c[i]; force=true) for i  in nm_c]
# [fix(b[i], data.b[i]; force=true) for i  in nm_b]
JuMP.register(m, :ev, 2*N, ev, autodiff=true)
@NLobjective(m, Max , sum(-ev(p[i,:] for i=1:N) for j=1:D)  )

unset_silent(m)
JuMP.optimize!(m)

st0 = termination_status(m)
obj0 = objective_value(m)
@show st0, obj0

psol0 = JuMP.value.(p)
csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)
zsol0 = [Asol0...;bsol0...;csol0...;psol0...]
Zsol0 = zsol0*zsol0'
rank(Zsol0)
eigvals(Zsol0)


cdf_pkg = [1-cdf(Beta(αsol0[i],βsol0[i]),w[i]/csol0[i]) for i=1:N]
cdf_opt = [dist_cdf([αsol0[i],βsol0[i],w[i]/csol0[i] ]...) for i=1:N]
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
    

    @test norm(cdf_pkg-cdf_opt)<0.01
    @test norm(delta-cdf_opt)<0.01
    
end

function clearing_vec!(F, p, x)
    F .= p.-min.(p_bar,max.((1+g0)*(Asol0'*p .+csol0 .- x.*csol0) .-g0.*p_bar,0))
end
function loss(x)
     # total loss in actual network for shock x
    p = nlsolve((F, p)->clearing_vec!(F, p, x),[p_bar...], autodiff = :forward)
    @test converged(p)
    sum(x.*csol0.+p_bar.-p.zero) 
end
xi = rand(dist,1000)
L_sim = mean(loss.(xi))
L_opt = sum(sum( (x0[i,j]*csol0[i]+p_bar[i]-psol0[i,j])*E.weights[j] for i=1:N) for j=1:D)
@test abs(L_sim-L_opt)/abs(L_opt)<0.01

Id = Matrix{Float64}(I,N,N)
p_lb = [min.(p_bar,(Id-Diagonal(x0[:,j]))*(Id-Asol0)\csol0) for j=1:D]
L_lb = sum(sum( (x0[i,j]*csol0[i]+p_bar[i]-p_lb[j][i])*E.weights[j] for i=1:N) for j=1:D)
L_d = sum(sum( (x0[i,j]*csol0[i]+max.(0,x0[i,j].-w[i]))*E.weights[j] for i=1:N) for j=1:D)

@show (L_opt,L_lb,L_d)

β = (p_bar.-bsol0)./p_bar # financial connectivity: proportion of liabilities inside the network
β⁺ = maximum(β)

ratio = L_opt/L_d
bound = 1 + sum(delta.*csol0)/((1-β⁺)*sum(csol0))
@show (ratio,bound)









# N=2
# dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))
# function dist_cdf(p...)
#     α,β,wc = p[1], p[2], p[3]
#     s(x,wc) = 1 ./(1 .+exp.(-1000*(x-wc))) # approx of indicator function of x>=wc
#     obj(x,p) = dist_pdf(x,p)*s(x,wc...)
#     prob = QuadratureProblem(obj,0.0,1.0,[α...,β...])
#     Quadrature.solve(prob,QuadGKJL(),reltol=1e-5,abstol=1e-5)[1]
# end
# p0 = [2.0,2.0,0.2]
# oo = dist_cdf(p0...)
# w = [3.0,3.0]
# delta = [0.03,0.23]
# m = Model(with_optimizer(Ipopt.Optimizer,max_iter=100000,print_level=5))
# JuMP.register(m, :dist_cdf, 3, dist_cdf,autodiff=true)
# @variable(m, 1.0<=α[i=1:N]<=3.0, start = 1.5) 
# @variable(m, 1.0<=β[i=1:N]<=50.0, start = 20.0) 
# @variable(m,c[i=1:N]>=w[i], start = 10.0*w[i]) 
# @NLexpression(m,wc[i=1:N],w[i]/c[i])
# for i in 1:N
#     @eval ($(Symbol("vv$i"))) = [α[$i],β[$i],wc[$i]]
#     @eval @NLconstraint(m,dist_cdf(($(Symbol("vv$i")))...)==delta[$i]) 
# end
# @NLobjective(m, Min, sum(wc[i] for i=1:N) )
# optimize!(m)
# a1 = JuMP.value.(α)
# b1 = JuMP.value.(β)
# c1 = JuMP.value.(c)
# termination_status(m)
# objective_value(m)
# cdf_pkg = [1-cdf(Beta(a1[i],b1[i]),w[i]/c1[i]) for i=1:N]
# cdf_opt = [dist_cdf([a1[i],b1[i],w[i]/c1[i] ]...) for i=1:N]
# @test norm(cdf_pkg-cdf_opt)<0.01
# @test norm(delta-cdf_opt)<0.01


using SparseGrids, FastGaussQuadrature

nodes, weights = gausslegendre( 100000 );
nodes2, weigths2 = sparsegrid(N, 1, gausslegendre)
g(x) = pdf(dist,x)
dot( weigths2, g(nodes2) )

using GalerkinSparseGrids
coeffs_DG(2,2, 5, x->pdf(dist,x); scheme="sparse")

dist = Product(Beta.([α0,β0]...)) 


Dirichlet(alpha) 

using ApproxFun

using PolyChaos

Beta01Measure
beta01
op = Uniform01OrthoPoly(n, addQuadrature=true);

degree, Nrec = 2, 20
α1, β1 = 1.0, 2.0
α2, β2 = 2.0, 10.0
op_beta1 = Beta01OrthoPoly(degree, α1, β1; Nrec=Nrec, addQuadrature=true)
op_beta2 = Beta01OrthoPoly(degree, α2, β2; Nrec=Nrec, addQuadrature=true)
mop=MultiOrthoPoly([op_beta1,op_beta2],degree)
L=PolyChaos.dim(mop)
f(t)=sin(t)
PolyChaos.integrate(f,mop)
