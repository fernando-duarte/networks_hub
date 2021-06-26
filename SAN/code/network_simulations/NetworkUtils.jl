module NetworkUtils

using Convex, SCS, NLsolve
using LinearAlgebra, Distributions, Random, SpecialFunctions,  SparseArrays, MeasureTheory
using BenchmarkTools, Test
using DataFrames, BSON, XLSX, JLD2, Missings
using Dates, MonthlyDates
# using Plots, Gadfly, Cairo, Fontconfig, LightGraphs, SimpleWeightedGraphs, GraphRecipes
using Setfield, Parameters, Suppressor, Flatten
using MathOptInterface
const MOI = MathOptInterface
using Reexport

include("IncBetaDer.jl")
include("NetworkType.jl")
@reexport using .NetworkType
using .IncBetaDer

export optimParam, c_lin, c_quad, c_chance, c_p, c_pNL, c_, c_NL, c_, c_NL, obj, obj_full, loss_d

struct optimParam{T}
    hᵀ::SparseMatrixCSC{T, Int64}
    q::Vector{T}
    H::SparseMatrixCSC{T, Int64}
    M::Integer
    net::Network 
    mini_batch::Integer
    selA::SparseMatrixCSC{T, Int64}
    selb::SparseMatrixCSC{T, Int64}
    selc::SparseMatrixCSC{T, Int64}
    selα::SparseMatrixCSC{T, Int64}
    selβ::SparseMatrixCSC{T, Int64}
    selp::SparseMatrixCSC{T, Int64}
    getA::Any
    getp::Any
    vars::Any
    #inner constructor    
    function optimParam(network::Network{T,N}, mini_batch::Integer=1) where {T,N}
        @unpack a, p_bar = network
        M = N^2+4*N+N*mini_batch
        #T = Base.typejoin(eltype(a),S)
        
        # constants for linear constraints 
        # hᵀz=q
        h₁₁ᵀ = hcat(sparse(kron(Matrix{T}(I,N,N),p_bar')))
        h₁ᵀ = hcat(h₁₁ᵀ,spzeros(T,N,N),sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
        q₁ = a

        h₂₁ᵀ = repeat(spdiagm(p_bar),1,N)
        h₂ᵀ = hcat(h₂₁ᵀ,sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
        q₂ = p_bar

        hᵀ = vcat(h₁ᵀ,h₂ᵀ)
        q = vcat(q₁,q₂)

        # constants for quadratic constraints
        # zᵀHz = 0
        li = LinearIndices((N,N))
        liᵀ= transpose(li)
        id = li[tril!(trues(N,N), 0)]
        idᵀ = liᵀ[tril!(trues(N,N), 0)]
        H = sparse(vcat(id,idᵀ),vcat(idᵀ,id),T(0.5),M,M)
        # add a linear term to the quadratic form
        # zᵀHz + h₃ᵀz = q₃
        #h₃ᵀ = spzeros(T,(1,M))
        #q₃ = spzeros(T,1)
        selA  = sparse(1:N^2,1:N^2,T(1),N^2,M)
        selb  = sparse(1:N,N^2+1:N^2+N,T(1),N,M)
        selc  = sparse(1:N,N^2+N+1:N^2+2N,T(1),N,M)
        selα  = sparse(1:N,N^2+2N+1:N^2+3N,T(1),N,M)
        selβ  = sparse(1:N,N^2+3N+1:N^2+4N,T(1),N,M)
        selp  = sparse(1:mini_batch*N,N^2+4N+1:M,T(1),mini_batch*N,M)
        # select variables in z
        getA(z) = reshape(selA*z,N,N)
        getp(z) = reshape(selp*z,N,mini_batch)
        vars(z) = (A=getA(z),b=selb*z,c=selc*z,α=selα*z,β=selβ*z,p=getp(z))

        new{T}(hᵀ,q,H,M,network,mini_batch,selA,selb,selc,selα,selβ,selp,getA,getp,vars)
    end
end #struct

## constraints and objective

function c_lin(z,op::optimParam)
    @unpack hᵀ, q  = op
    return hᵀ*z - q  #hᵀz=q
end 

function c_quad(z,op::optimParam)
    @unpack H  = op
    return dot(z,H*z) #transpose(z)*H*z
end

## chance constraints ⟺ c_chance(z)=0
# Prob(c[i]x[i]>=w[i])=delta[i] 
beta_ccdf(a,b,x) = zero(x) < x < one(x) ? one(x) - beta_inc_grad(a, b, x)[1] : zero(x)
function c_chance(z,op::optimParam)
    @unpack selα, selβ, selc, net  = op
    @unpack w, delta = net
    zeroz = zero(z[1])+0.0001
    onez = one(z[1])-0.0001
    return beta_ccdf.(selα*z,selβ*z,max.(min.(w./(selc*z),onez),zeroz) ) .- delta
end

## clearing vector constraint
# quadratic version (ignores min{p_bar,max{0,...}} )
function c_p(z,x,op::optimParam)
    @unpack selc, getA, getp, net  = op 
    @unpack p_bar, γ, N = net
    c = selc*z
    p = getp(z)
    Aᵀ = transpose(getA(z))
    return reshape(p - ((1+γ)*Aᵀ*p + ((1+γ)*(one(x[1]) .- x)).*repeat(c,1,size(x,2))-γ*repeat(p_bar,1,size(x,2))) ,N*size(x,2))
end

# non-linear version
function c_pNL(z,x,op::optimParam)
    @unpack selc, getA, getp, net  = op 
    @unpack p_bar, γ, N = net
    c = selc*z
    p = getp(z)
    Aᵀ = transpose(getA(z))
    p_bar_rep = repeat(p_bar,1,size(x,2))
    return reshape(p-min.(max.(0,(1+γ)*Aᵀ*p + ((1+γ)*(one(x[1]) .- x)).*repeat(c,1,size(x,2))-γ*p_bar_rep),p_bar_rep),N*size(x,2))
end

# write constraints as cL<=c(z)<=cU, setting cL=cU for equality constraints
function c_(z,x,op::optimParam)
    #x = rand(N,mini_batch)
    #x = ones(N,mini_batch)/10
    vcat(c_lin(z,op),
        c_quad(z,op),
        c_chance(z,op),
        c_p(z,x,op)
        )
end

function c_NL(z,x,op::optimParam)
    #x = rand(N,mini_batch)
    #x = ones(N,mini_batch)/10
    vcat(c_lin(z,op),
        c_quad(z,op),
        c_chance(z,op),
        c_pNL(z,x,op)
        )
end

function c_(z,op::optimParam)
    vcat(c_lin(z,op),
        c_quad(z,op),
        c_chance(z,op)
        )
end

function c_NL(z,op::optimParam)
    vcat(c_lin(z,op),
        c_quad(z,op),
        c_chance(z,op)
        )
end

# objective to max or min
function obj(z,op::optimParam)
    @unpack selc, selα, selβ, getp  = op 
    c = selc*z
    α = selα*z
    β = selβ*z
    p = getp(z)
    return sum(c.*α./(α + β))-sum(mean(p,dims=2))
end

function obj(z,x,op::optimParam)
    @unpack selc, selα, selβ, getA, getp, net = op 
    @unpack p_bar, γ = net
    Aᵀ = transpose(getA(z))
    c = selc*z
    α = selα*z
    β = selβ*z
    mini_batch = size(x,2)
    p_guess = getp(z) 
    p_guess_rep = repeat(p_guess,1,mini_batch)
    p_bar_rep = repeat(p_bar,1,mini_batch)
    c_rep = repeat(c,1,mini_batch)
    function clearing_p!(F, p, x, Aᵀ, c,γ,p_bar_rep,c_rep)
        F .= p.-min.(p_bar_rep,max.( (1+γ)*(Aᵀ*p .+ (one(x[1]) .- x).*crep) .- γ.*p_bar_rep, 0 ) )
    end
    p = nlsolve((F, p)->clearing_p!(F, p, x, Aᵀ, c,γ,p_bar_rep,c_rep,p_guess_rep), autodiff = :forward)
    #@test converged(p)
    #return p.zero
    return sum(c.*x)-sum(mean(p.zero,dims=2)) # spillovers for network 
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end

function obj_full(z,op::optimParam)
    @unpack net  = op
    @unpack p_bar = net
    obj(z,op)+sum(p_bar)
end

function obj_full(z,x,op::optimParam)
    @unpack net  = op
    @unpack p_bar = net
    obj(z,x,op)+sum(p_bar)
end

function loss_d(z,x,op::optimParam)
    @unpack selc, net = op
    @unpack w = net
    c = selc*z
    cx = c.*x
    return sum(cx+max.(0,cx-w)) # spillovers for disconnected network 
end



end #module