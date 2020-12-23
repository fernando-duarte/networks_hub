## load packages
using Pkg
Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("NLsolve")
Pkg.add("Complementarity")
Pkg.add("NLopt")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("SpecialFunctions")
Pkg.add("ForwardDiff")
Pkg.add("AmplNLWriter")
Pkg.add("Cbc")

Pkg.add("FFTW")
Pkg.build("FFTW")
ENV["GRDIR"]="" ; Pkg.build("GR")
Pkg.add("PATHSolver")
Pkg.build("PATHSolver")
Pkg.add("SLEEFPirates")
Pkg.build("SLEEFPirates")
Pkg.add("IterTools")
Pkg.add("Tracker")
Pkg.add("Zygote")
Pkg.add("SpecialFunctions")

Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("NLsolve")
Pkg.add("Complementarity")
Pkg.add("NLopt")
Pkg.add("JuMP")
Pkg.add("Ipopt")

Pkg.add("ForwardDiff")
Pkg.add("GraphRecipes")
Pkg.add("Gadfly")

Pkg.add("LightGraphs")
Pkg.add("SimpleWeightedGraphs")
Pkg.add("CUDA")
Pkg.add("CuArrays")
Pkg.add("Distributed")
Pkg.add("Surrogates")
Pkg.add("QuadGK")
Pkg.add("DistributionsAD")
Pkg.add("StatsFuns")

Pkg.rm("CuArrays")
Pkg.rm("DiffEqFlux")
Pkg.rm("NeuralNetDiffEq")
Pkg.add(Pkg.PackageSpec(url="https://github.com/FluxML/Flux.jl", rev="master"))
Pkg.update("CUDA")
Pkg.build("CUDA")
Pkg.add("DistributionsAD")

Pkg.add(Pkg.PackageSpec(url="https://github.com/SciML/DiffEqFlux.jl", rev="master"))


using LinearAlgebra, DataFrames, XLSX, Missings
using Distributions, DistributionsAD
using Random
using CUDA
using Flux
using Flux.CUDA
using Flux: cpu, gpu
using Flux: @epochs
using Flux.Data 
using Test
using IterTools: ncycle
using Zygote

# using GraphRecipes
# using Plots
# using Gadfly
# using LightGraphs
# using SimpleWeightedGraphs

# CUDA.allowscalar(false)



#= network primitives from the data
p_bar = total liabilities
c = outside assets
assets = total assets
w = net worth
b = outside liabilities
=#

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... )
unique!(data) # delete duplicate rows, use `nonunique(data)' to see if there are any duplicates
#data[(isequal.(data.tkr,"JPM") .| isequal.(data.tkr,"AIG")) .& isequal.(data.qt_dt,192),:]
#data = data[isless.(0,data.w), :]
data = data[isequal.(data.qt_dt,192), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:50,:] # keep small number of nodes, for testing
N = size(data,1) # number of nodes

# rescale units
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units

# fake missing data to make problem feasible
data.b[:] .= missing
data.c[:] .= missing

# create network variables
#data.f = data.p_bar .- data.b # inside liabilities
#data.d =  data.assets .- data.c # inside assets

# keep track of missing variables
col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing
#data_nm = dropmissing(data, disallowmissing=true) # drop missing
data_nm = coalesce.(data, data.w .+ 0.01) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))
dropmissing(data, [:delta, :delta_alt, :w, :assets, :p_bar]) # remove type missing

# take a look
names(data) # column names
describe(data)
show(data, true)

# function shock(c)
#     # draw shocks that produce the probability of default δ
#     a = 1
#     b = log.(data.delta)./(log.(1 .-data.w./c))
#     dist = Beta.(a,[b...])
#     draws = rand.(dist, 1)
#     vcat(draws...).*c
#     #rand(Float64,N) # normal
# end


## Optimization
# parameters
CUDA.allowscalar(true)

temp = convert(Array{Float32},data[:,[:delta, :delta_alt, :w, :assets, :p_bar]])
delta = copy(temp[:,1]); delta = gpu(delta)
delta_alt = copy(temp[:,2]); delta_alt = gpu(delta_alt)
w= copy(temp[:,3]); w = gpu(w)
assets= copy(temp[:,4]); assets = gpu(assets)
p_bar = copy(temp[:,5]); p_bar = gpu(p_bar)

g0 = Float32[0.0] |> gpu # bankruptcy cost
rng =  Random.seed!(123)
α=1.0f0  |>gpu;β=2.0f0  |>gpu;
dist = Beta(α,β)
dist = convert(Beta{Float32},dist) |> gpu

D= 1
x =  convert(Array{Float32},rand.(dist, N, D)) |> gpu
σA = convert(Array{Float32},rand(rng,N,N))  |> gpu
σA .= -0.9
σA[diagind(σA)] .= -1 
for i=1:N
    for j=i+1:N
        σA[i,j] = -1 
    end
end

σA .= -10
σA[diagind(σA)] .= -10
for i=1:N
    for j=i+1:N
        σA[i,j] = -10
    end
end
σc = fill(1.0f0,N) |>gpu
σb = fill(0.0f0,N) |>gpu
# σA = σA |> gpu


## fixed point
function contraction(p,x)
    #A = (NNlib.hardtanh.(σA).+cu(1))./cu(2)
    A = NNlib.sigmoid.(σA)
    c = NNlib.sigmoid.(σc).*assets
    #c = -A' * p_bar .+ assets 
   min.(p_bar, max.((1 .+g0).*(A'*p .+ c .- x.*c) .- g0.*p_bar,0))
end

# need to modify for x with more than one column
contraction_iter(x, n::Integer) = n <= 0 ? p_bar  : contraction(contraction_iter(x,n-1),x)

function betapdf(x) 
    #A = (NNlib.hardtanh.(σA).+cu(1))./cu(2)
    A = NNlib.sigmoid.(σA)
    Aᵀ = NNlib.sigmoid.(transpose(σA))
    c = NNlib.sigmoid.(σc).*assets
    #c = -A' * p_bar .+ assets 
    #α = 1.0f0 |> gpu
    wc = min.(max.(1f-6, w./c),1-1f-6)
    β=log.(delta)./(log.(1 .- wc))
    #x.^(α .- 1) .* (1 .- x).^(β .- 1) 
    (1 .- x).^(β .- 1) 
end

loss_x(x) = -sum(contraction_iter(x,3).*betapdf(x),dims=1)
#loss_x(x) = -sum(contraction_iter(x,3),dims=1)

# Δx = 0.01
# gridx = 0:Δx:1
# function loss_int(σA,σb,σc)
#     A = (NNlib.hardtanh.(σA).+1)./2
#     b = NNlib.sigmoid.(σb).*p_bar
#     c = NNlib.sigmoid.(σc).*assets
#     sum(loss_x.(gridx,Ref(σA),Ref(σb),Ref(σc)).*exp.(logpdf.(dist,gridx))*Δx )
# end
cu_ones = ones(Float32,N) |> gpu
function penalty(x)
#     A = (NNlib.hardtanh.(σA).+1)./2
#     Aᵀ= (NNlib.hardtanh.(transpose(σA)).+1)./2
    A = NNlib.sigmoid.(σA)
    Aᵀ = NNlib.sigmoid.(transpose(σA))
    b = NNlib.sigmoid.(σb).*p_bar
    c = NNlib.sigmoid.(σc).*assets
    con1 = sum(A,dims=2).*p_bar - p_bar .+ b
    con2 = A' * p_bar .- assets .+ c
    con3 = sum(LowerTriangular(A).*LowerTriangular(Aᵀ)) 
    #sum(LowerTriangular(A).*LowerTriangular(Aᵀ))==sum(LowerTriangular(A.*Aᵀ)) check which is faster
     sum(abs2,con1)+sum(abs2,con2)+con3.^2
    #con3.^2
end

min_obj = 1 # minimize (min_obj=1) or maximize (min_obj=-1)
function lossFlux(x)
   min_obj*mean(loss_x(x)) + penalty(x)
end

ps = Flux.params(σA,σc) 
pabstol = 1e-8

cbP = function ()
    #l < pabstol && Flux.stop() 
    display((lossFlux(x),min_obj*mean(loss_x(x)) , penalty(x)))
end
cbPenalty = function ()
    display(penalty(x))
end
cbPenalty()
cbP()

train_loader = DataLoader(x, batchsize=1, shuffle=true) 
# opt = Flux.ADAM(0.1)
# Flux.train!(penalty, ps ,ncycle(train_loader, 30) , opt, cb = cbPenalty)
opt = Flux.ADAM(0.01)
@epochs 10 Flux.train!(lossFlux, ps ,ncycle(train_loader, 1) , opt, cb = cbP)

function my_custom_train!(loss, ps, data , opt)
  # training_loss is declared local so it will be available for logging outside the gradient calculation.
  local training_loss
  ps = Params(ps)
  data = DataLoader(CUDA.rand(N,10), batchsize=1, shuffle=true)
  train_data = ncycle(data, 1)
  for d in train_data
    gs = gradient(ps) do
      training_loss = loss(d)
      # Code inserted here will be differentiated, unless you need that gradient information
      # it is better to do the work outside this block.
      return training_loss
    end
    # Insert whatever code you want here that needs training_loss, e.g. logging.
    # logging_callback(training_loss)
    # Insert what ever code you want here that needs gradient.
    # E.g. logging with TensorBoardLogger.jl as histogram so you can see if it is becoming huge.
    Flux.update!(opt, ps, gs)
    # Here you might like to check validation set accuracy, and break out to do early stopping.
    display((lossFlux(x),min_obj*mean(loss_x(x)) , penalty(x)))
  end
end
@epochs 1000 my_custom_train!(lossFlux, ps ,ncycle(train_loader, 1) , opt)
@epochs 300 my_custom_train!(penalty, ps ,ncycle(train_loader, 1) , opt)


for d in 
    #display(d)
    display(lossFlux(d))
end

#Asol = (NNlib.hardtanh.(σA).+1)./2
Asol = NNlib.sigmoid.(σA) |> cpu
csol = NNlib.sigmoid.(σc).*assets |> cpu
bsol = NNlib.sigmoid.(σb).*p_bar |> cpu

tol = 1e-4
@testset "check solution" begin
   @test norm( sum(Asol,dims=2).* p_bar .- (p_bar .- bsol)) < tol
   @test norm( Asol' * p_bar .- (assets .- csol)) < tol
   @test norm(diag(Asol)) < tol
   @test norm([Asol[i,j]*Asol[j,i] for i=1:N , j=1:N]) < tol
   @test all(0 .<=Asol.<=1)
   @test all(0 .<=bsol.<=p_bar)
   @test all(0 .<=csol.<=assets)
end
Aplot = copy(Asol)
Aplot[Aplot.<1f-4].=0

Gadfly.spy(Asol)

function clearing_vec!(F, p, x, A)
    F .= p.-min.(p_bar,max.((1+g0)*(A'*p .+c .- x) .-g0.*p_bar,0))
end

function loss_c(x,A)
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

p_bar = p_bar |> cpu
delta = delta |> cpu
β = β |> cpu

β = (p_bar.-bsol)./p_bar # financial connectivity: proportion of liabilities inside the network
β⁺ = maximum(β)
ratio_x = loss_c(x,Asol)/loss_d(x)
ratio = mean(loss_x(x))
bound = 1 + sum(delta.*csol)/((1-β⁺)*sum(csol))


# Pkg.add("JLD")
# using JLD
# save("/home/ec2-user/SageMaker/Test-AWS/net_opt.jld", "Asol", Asol,"data",data)

## plot
Aplot = deepcopy(Asol)
Aplot[Aplot.<1e-3] .=0

# attributes here: https://docs.juliaplots.org/latest/generated/graph_attributes/
#method `:spectral`, `:sfdp`, `:circular`, `:shell`, `:stress`, `:spring`, `:tree`, `:buchheim`, `:arcdiagram` or `:chorddiagram`.


graphplot(LightGraphs.DiGraph(Aplot),
         nodeshape=:circle,
         markersize = 0.05,
         node_weights = assets,
         markercolor = range(colorant"yellow", stop=colorant"red", length=N),
         names = data.nm_short,
         fontsize = 8,
         linecolor = :darkgrey,
         edgewidth = (s,d,w)->5*Asol[s,d],
         arrow=true,
         method= :circular, #:chorddiagram,:circular,:shell
         )


## some tests
tol = 1e-6
x_test = fill(0.0,N,1)
pinf = nlsolve(p->contraction(p, x_test)-p, p_bar, autodiff = :forward)
@test norm(contraction_iter(x_test,100)-pinf.zero)<tol

Zygote.gradient(c -> -sum(contraction(p_bar,x_test)),data_nm.c)
ForwardDiff.gradient(c -> -sum(contraction(p_bar,x_test)),data_nm.c)

Zygote.gradient(c -> -sum(contraction_iter(x_test,20)),data_nm.c)
ForwardDiff.gradient(c -> -sum(contraction_iter(x_test,20)),data_nm.c)

Zygote.gradient(loss_x(x_test),data_nm.c)
Zygote.gradient(c -> -sum(contraction_iter(x0,3,c).*exp.(logpdf.(dist,0.0))),data_nm.c)

loss(data_nm.c)
zy=Zygote.gradient(loss, data_nm.c)
fd=ForwardDiff.gradient(loss, data_nm.c)



## debug code

for d in train_data
    display(d)
end

Zygote.gradient(()->penalty(x),σb)

ps = Params(ps)
  for d in xx
    gs = Zygote.gradient(ps) do
      training_loss = penalty(d...)
      # Code inserted here will be differentiated, unless you need that gradient information
      # it is better to do the work outside this block.
      return training_loss
    end
end
    # Insert whatever code you want here that needs training_loss, e.g. logging.
    # logging_callback(training_loss)
    # Insert what ever code you want here that needs gradient.
    # E.g. logging with TensorBoardLogger.jl as histogram so you can see if it is becoming huge.
    #update!(opt, ps, gs)
    # Here you might like to check validation set accuracy, and break out to do early stopping.
  #end





# x1 = Float32[0.024476904, 0.19392614, 0.4850096] |> gpu
# x2 = Float32[0.09934894, 0.09108944, 0.3968748] |> gpu
# train_loader = [(x1,),(x2,)]
# @epochs 3 Flux.train!(loss_x, ps ,train_loader , opt, cb = cbP)

# train_loader = DataLoader((x,), batchsize=1) 
# @epochs 3 Flux.train!(loss_x, ps ,train_loader , opt, cb = cbP)
# Flux.train!(loss_x, ps ,ncycle(train_loader, 3) , opt, cb = cbP)
# @epochs 3 Flux.train!(loss_x, ps ,train_loader , opt, cb = cbP)

import Pkg; Pkg.add("Complementarity")
using JuMP
using Ipopt
using Complementarity 
using Flux.NNlib

A0 = (NNlib.hardtanh.(σA).+1)./2 |>cpu
b0 = NNlib.sigmoid.(σb).*data.p_bar |>cpu
c0 = NNlib.sigmoid.(σc).*data.assets |>cpu

# set up optimization
m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps",max_iter=10000))

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = c0[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = b0[i])   
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  # start=A0[i,j]

@constraint(m, sum(A,dims=2).*data.p_bar .==  data.p_bar .- b) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c) # payments from other nodes add up to inside assets d

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
    push!(myexpr, (1+g0[1])*(A'*p[:,j] .+ c .- x[:,j]) .- g0[1]*data.p_bar)
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

#@NLobjective(m, Max , sum(sum( x0[i,j]*c[i]+data.p_bar[i]-p[i,j] for i=1:N)/D for j=1:D) ) #*sum(x[i]+p_bar[i]-p[i] for i=1:N) 
@NLobjective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 

#unset_silent(m)
JuMP.optimize!(m)

termination_status(m)
objective_value(m)

psol = JuMP.value.(p)
csol = JuMP.value.(c)
bsol = JuMP.value.(b)
Asol = JuMP.value.(A)

tol = 1e-6
@testset "check solution" begin
   @test norm( sum(Asol,dims=2).* p_bar .- (p_bar .- bsol)) < tol
   @test norm( Asol' * p_bar .- (assets .- csol)) < tol
   @test norm(diag(Asol)) < tol
   @test norm([Asol[i,j]*Asol[j,i] for i=1:N , j=1:N]) < tol
   @test all(0 .<=Asol.<=1)
   @test all(0 .<=bsol.<=p_bar)
   @test all(0 .<=csol.<=assets)
end

Asol[Asol.<1f-4].=0


#import Pkg; Pkg.add("Optim"); Pkg.add("NLSolversBase")
using Optim, NLSolversBase
import NLSolversBase: clear!
using Zygote

x = x |> cpu
p_bar = p_bar |> cpu
g0 = g0 |> cpu
assets = assets |> cpu
b0 = -sum(A0,dims=2).*p_bar .+ p_bar 
c0 = -A0' * p_bar .+ assets 
A0[1,1]=0.1 
A0 = rand(N,N)
θ₀ = vcat(A0...,b0,c0) 
θ₀ = [θ₀...]
θ₀ = convert.(Float32,θ₀)
 
x = convert.(Float64,x)
p_bar = convert.(Float64,p_bar)
g0 = convert.(Float64,g0)
assets = convert.(Float64,assets)
θ₀ = convert.(Float64,θ₀)


function contraction(p,x,θ)
    A = reshape(θ[1:N^2],N,N)
    b = θ[N^2+1:N^2+N]
    c = θ[N^2+N+1:end]
    min.(p_bar, max.((1 .+g0[1]).*(A'*p .+ c .- x.*c) .- g0[1].*p_bar,0))
end

contraction_iter(x, n::Integer,θ) = n <= 0 ? p_bar  : contraction(contraction_iter(x,n-1,θ),x,θ)
α=1.0  ;β=2.0 ;
betapdf(x) = x^(α - 1) * (1 - x)^(β - 1) 
loss_x(x,θ) = mean(sum(contraction_iter(x,1,θ).*betapdf.(x),dims=1))

function penalty(θ)
    A = reshape(θ[1:N^2],N,N)
    b = θ[N^2+1:N^2+N]
    c = θ[N^2+N+1:end]
    con1 = sum(A,dims=2).*p_bar .- p_bar .+ b
    con2 = A' * p_bar .- assets .+ c
    con3 = sum(LowerTriangular(A).*LowerTriangular(transpose(A)))
    vcat(con1...,con2...,con3...)
end

fun(θ) = loss_x(x,θ)
gs = Zygote.gradient(Params(θ₀)) do
         fun(θ₀)
end
gs[θ₀]

function fun_grad!(g, θ)
    gtemp = Zygote.gradient(fun,θ)
    g=gtemp[1]
    #g
end   
function fun_grad(θ)
    g = Zygote.gradient(fun,θ)[1]
    g
end   

function fun_hess!(h, θ)
    h = Zygote.hessian(fun,θ)
    #h
end
function fun_hess(θ)
    h = Zygote.hessian(fun,θ)
    h
end

function con_c!(c, θ)
    c = penalty(θ)
    #c
end
function con_c(θ)
    c = penalty(θ)
    c
end

function con_jacobian!(J, θ)
    A = reshape(θ[1:N^2],N,N)
    dcon1 = hcat(transpose(kron(Matrix{Float64}(I,N,N),p_bar)),Matrix{Float64}(I,N,N),fill(0.0 ,N,N))
    dcon2 = hcat(transpose(kron(Matrix{Float64}(I,N,N),p_bar)),fill(0.0 ,N,N),Matrix{Float64}(I,N,N))
    dcon3 = hcat((transpose(A)+Diagonal(A))...,fill(0.0 ,1,N),fill(0.0 ,1,N)) 
    J = vcat(dcon1,dcon2,dcon3)
    #J
end  

function con_jacobian(θ)
    A = reshape(θ[1:N^2],N,N)
    dcon1 = hcat(transpose(kron(Matrix{Float64}(I,N,N),p_bar)),Matrix{Float64}(I,N,N),fill(0.0 ,N,N))
    dcon2 = hcat(transpose(kron(Matrix{Float64}(I,N,N),p_bar)),fill(0.0 ,N,N),Matrix{Float64}(I,N,N))
    dcon3 = hcat((transpose(A)+Diagonal(A))...,fill(0.0 ,1,N),fill(0.0 ,1,N)) 
    J = vcat(dcon1,dcon2,dcon3)
    J
end 

Hconst = zeros(N^2,N^2);
for i=1:N
    for j=1:N
        temp = zeros(N,N)
        if i==j
            temp[j,i]=2.0 
        else
            temp[j,i]=1.0 
        end
        Hconst[(i-1)*N+1:i*N,(j-1)*N+1:j*N] = temp
    end
end
function con_h!(h, θ, λ)
    # con3
    h[1:N^2,1:N^2] = h[1:N^2,1:N^2]+ λ[2*N+1]*Hconst
    #h
end
function con_h(θ, λ)
    # con3
    h[1:N^2,1:N^2] = h[1:N^2,1:N^2]+ λ[2*N+1]*Hconst
    h
end

fun(θ₀)
typeof(fun_grad!(fill(1.0 ,length(θ₀)),θ₀))
typeof(fun_hess!(fill(1.0 ,length(θ₀),length(θ₀)),θ₀))
con_c!(fill(1.0 ,2*N+1),θ₀)
con_jacobian!(fill(1.0 ,length(θ₀)),θ₀)
con_h!(fill(1.0 ,length(θ₀),length(θ₀)),θ₀,fill(1.0 ,2*N+1))

clear!(df)
df =TwiceDifferentiable(fun, fun_grad!, fun_hess!, rand(length(θ₀)))
lx = fill(0.0 ,length(θ₀)); ux = vcat(fill(1.0 ,N^2),p_bar,assets)
lc = fill(-1.0 ,2*N+1); uc = fill(1.0 ,2*N+1)
dfc = TwiceDifferentiableConstraints(con_c!, con_jacobian!, con_h!,lx, ux, lc, uc)
res = optimize(df, dfc,rand(length(θ₀)), IPNewton())


clear!(df)
lx = fill(0.0 ,length(θ₀)); ux = vcat(fill(1.0 ,N^2),p_bar,assets)
c = fill(-1000 ,2*N+1); uc = fill(1000 ,2*N+1)
t0 = rand(length(θ₀)).*(ux/2)
df =TwiceDifferentiable(fun, fun_grad, fun_hess,t0, inplace = false)
dfc = TwiceDifferentiableConstraints(con_c!, con_jacobian!, con_h!,lx, ux, lc, uc)
#dfc = TwiceDifferentiableConstraints(lx, ux)
res = optimize(df, dfc,t0, IPNewton())

res = optimize(df, dfc,θ₁, IPNewton())

θ₁ = Optim.minimizer(res)
Optim.minimum(res)

df =TwiceDifferentiable(fun, fun_grad, fun_hess, θ₀, inplace = false)
dfc = TwiceDifferentiableConstraints(lx, ux)
optimize(df, dfc, θ₀, IPNewton())
