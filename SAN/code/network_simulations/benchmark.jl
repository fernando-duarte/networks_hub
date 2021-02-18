import Pkg; Pkg.add("Pkg")
using Pkg
Pkg.build("HDF5")
Pkg.add("JLD")
Pkg.add("SpecialFunctions")
Pkg.build("SpecialFunctions")
Pkg.add("FFTW")
Pkg.build("FFTW")
Pkg.build("PATHSolver")
Pkg.add("PATHSolver")
Pkg.build("SLEEFPirates")
Pkg.add("SLEEFPirates")
Pkg.add("IterTools")
Pkg.add("Tracker")
Pkg.build("DiffEqBase")
Pkg.add("DiffEqBase")

# restart kernel after running the lines above
using Pkg
Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("Complementarity")
Pkg.add("DistributionsAD")
Pkg.add("Zygote")
Pkg.add(Pkg.PackageSpec(url="https://github.com/SciML/Quadrature.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaDiff/ForwardDiff.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaDiff/FiniteDiff.jl", rev="master"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/giordano/Cuba.jl", rev="master"))
Pkg.add("Cubature")
Pkg.add("HCubature")
Pkg.add("PlotlyBase");Pkg.add("PlotlyJS"); Pkg.add("ORCA")
 
using LinearAlgebra, DataFrames, XLSX, Missings, JuMP, Ipopt, Random, Complementarity, Test, Distributions, DistributionsAD, SpecialFunctions, NLsolve
using Quadrature, ForwardDiff, FiniteDiff, Zygote, Cuba, Cubature, HCubature
using Plots, Profile
plotly(ticks=:native)  

N = 20 # keep largest `N' nodes by assets

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
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
β0= []
for i=1:N
    function ff!(F, x)
        F[1] = 1-cdf(Beta(α0[i],x[1]),w[i]/c0[i])-delta[i] 
    end
    if i==1
        sol = nlsolve(ff!,[50.0])
    else
        sol = nlsolve(ff!,[1.0])
    end
    @test converged(sol)
    push!(β0,sol.zero[1])
end

dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))
Plots.plot([x->dist_pdf(x,[α0[i],β0[i]]) for i=1:N],0,1)

function dist_cdf(p...)
    α,β,wc = p[1], p[2], p[3]
    s(x,wc) = 1 ./(1 .+exp.(-1000*(x-wc))) # approximates indicator function of x>=wc
    obj(x,p) = dist_pdf(x,p)*s(x,wc...)
    prob = QuadratureProblem(obj,0.0,1.0,[α...,β...])
    Quadrature.solve(prob,QuadGKJL(),reltol=1e-5,abstol=1e-5)[1]
end

function dist_pdf_mv(x,params)
           # α,β,A,b,c = params[1][1:N], params[1][N+1:2*N], params[1][2*N+1:2*N+N^2], params[1][2*N+N^2+1:3*N+N^2], params[1][3*N+N^2+1:4*N+N^2]
            α,β,A,b,c = params[1:N], params[N+1:2*N], params[2*N+1:2*N+N^2], params[2*N+N^2+1:3*N+N^2], params[3*N+N^2+1:4*N+N^2]
            A = reshape([A...],N,N) 
            joint_pdf = 1.0;
            for i=1:N
               joint_pdf = joint_pdf*dist_pdf(x[i],[α[i],β[i]]) 
            end

            function contraction(x,p)
               min.(p_bar, max.((1 .+g0).*(A'*p .+ c .- x.*c) .- g0.*p_bar,0))
            end
            contraction_iter(x, n::Integer) = n <= 0 ? p_bar  : contraction(x,contraction_iter(x,n-1))
            loss_x(x) = -sum(contraction_iter(x,2)).*dist_pdf(x,[α...,β...])   
            loss_x(x)*joint_pdf
end


# function dist_pdf_mv_test(x,params...)
#     α,β,A,b,c = params[1:N], params[N+1:2*N], params[2*N+1:2*N+N^2], params[2*N+N^2+1:3*N+N^2], params[3*N+N^2+1:4*N+N^2]
#     A = reshape([A...],N,N) 
#     joint_pdf=1.0;
#     for i=1:N
#         joint_pdf = joint_pdf*dist_pdf(x[i],[α[i],β[i]]) 
#     end
#     #joint_pdf = prod(dist_pdf.(x,[α,β]))
#     some_fun = x->sum(x) # some arbitrary function to test
#     some_fun([α...,β...,A...,b...,c...])*joint_pdf 
# end
# dist_pdf_mv_test(ones(N)/12.0,ic...)

function ev(params...)
    prob = QuadratureProblem(dist_pdf_mv,zeros(N),ones(N),[params...])
    Quadrature.solve(prob,CubaCuhre(),reltol=1e-5,abstol=1e-5)[1]
end

ic = [α0...,β0...,A0...,b0...,c0...]
dist_pdf_mv(ones(N)/10.0,ic)
ev(ic...)

# set up optimization
# m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps",max_iter=10000,hessian_approximation="limited-memory"))
m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>100000,"print_level"=>5))

@variable(m, 0.8<=α[i=1:N]<=3.0, start = α0[i]) 
@variable(m, 1.0<=β[i=1:N]<=150.0, start = β0[i]) 
#@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, w[i]<=c[i=1:N]<=data.assets[i], start = c0[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = b0[i])   

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

# match probabilities of default
JuMP.register(m, :dist_cdf, 3, dist_cdf, autodiff=true)
@NLexpression(m,wc[i=1:N],w[i]/c[i])
for i in 1:N
    @eval ($(Symbol("vv$i"))) = [α[$i],β[$i],wc[$i]]
   @eval @NLconstraint(m,dist_cdf(($(Symbol("vv$i")))...)== delta[$i]) 
end

#[fix(c[i], data.c[i]; force=true) for i  in nm_c] #fixing c to data
#[fix(b[i], data.b[i]; force=true) for i  in nm_b] #fixing b to data
JuMP.register(m, :ev, 4*N+N^2, ev, autodiff=true)
vars = [α...,β...,A...,b...,c...]
@NLobjective(m, Max , ev(vars...)  )

unset_silent(m)
@time JuMP.optimize!(m)

st0 = termination_status(m)
obj0 = objective_value(m)
@show st0, obj0

csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)

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

    @test norm(cdf_pkg-cdf_opt)<0.01
    @test norm(delta-cdf_opt)<0.01
    
end

# Displaying output in batch solution




#=
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

=#







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

#=


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





dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))

function dist_pdf_mv2(x,params)
     α,β,A,b,c = params[1:N], params[N+1:2*N], params[2*N+1:2*N+N^2], params[2*N+N^2+1:3*N+N^2], params[3*N+N^2+1:4*N+N^2]
    A = reshape([A...],N,N) 
    joint_pdf = 1.0;
    for i=1:N
       joint_pdf = joint_pdf*dist_pdf(x[i],[α[i],β[i]]) 
    end
    
    function contraction(x,p)
       min.(p_bar, max.((1 .+g0).*(A'*p .+ c .- x.*c) .- g0.*p_bar,0))
    end
    contraction_iter(x, n::Integer) = n <= 0 ? p_bar  : contraction(x,contraction_iter(x,n-1))
    loss_x(x) = -sum(contraction_iter(x,2)).*dist_pdf(x,[α...,β...])   
    loss_x(x)*joint_pdf
end

function ev2(params)
    prob2 = QuadratureProblem(dist_pdf_mv2,zeros(N),ones(N),params)
    Quadrature.solve(prob2,CubaCuhre(),reltol=1e-5,abstol=1e-5)[1]
end


dist_pdf_mv(ones(N)/5,ic)
dist_pdf_mv2(ones(N)/5,ic)
ev(ic...)
ev2(ic)

ForwardDiff.gradient(p->dist_pdf_mv(ones(N)/3,p), ic )
ForwardDiff.gradient(p->dist_pdf_mv2(ones(N)/3,p), ic )
ForwardDiff.gradient(p->ev(p), ic )
∇ev2(params) = ForwardDiff.gradient(p->ev2(p), params )
∇ev2(ic)

ic = [α0...,β0...,A0...,b0...,c0...]
dist_pdf_mv(ones(N)/20,[ic...])
ev(ic...)
ev(ic)
dev = ForwardDiff.gradient(dist_pdf_mv,params)


=#