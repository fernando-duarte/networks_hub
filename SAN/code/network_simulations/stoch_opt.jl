using LinearAlgebra, DataFrames, XLSX, Missings, Test, SpecialFunctions,  SparseArrays , Random, BenchmarkTools, Distributions
using JuMP,Ipopt, NLsolve, Optim, LineSearches
using NLPModels, ADNLPModels, Percival, Nonconvex, Hyperopt
using ForwardDiff
using InteractiveErrors
include("IncBetaDer.jl")

## load data

xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
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

N=5
T = Float64
temp = Array{T}(data[1:N,[:delta, :delta_alt, :w, :assets, :p_bar]])
delta = copy(temp[:,1]); 
delta_alt = copy(temp[:,2]);
w = copy(temp[:,3]); 
assets= copy(temp[:,4]);
p_bar = copy(temp[:,5]);

γ = T(0)
mini_batch = 10;
M = N^2+4*N+mini_batch*N

# z = vcat(A...,b,c,α,β,p...)
# p is N-by-mini_batch
@noinline function main()
    selA  = sparse(1:N^2,1:N^2,T(1),N^2,M)
    selb  = sparse(1:N,N^2+1:N^2+N,T(1),N,M)
    selc  = sparse(1:N,N^2+N+1:N^2+2N,T(1),N,M)
    selα  = sparse(1:N,N^2+2N+1:N^2+3N,T(1),N,M)
    selβ  = sparse(1:N,N^2+3N+1:N^2+4N,T(1),N,M)
    selp  = sparse(1:mini_batch*N,N^2+4N+1:M,T(1),mini_batch*N,M)

    getA(z) = reshape(selA*z,N,N)
    getp(z) = reshape(selp*z,N,mini_batch)

    ## upper and lower bounds
    # low<=z<=upp
    uppα = 10
    uppβ = 50

    low_zero = eps(T) # eps(T) T(0.005)
    low = vcat(spzeros(T,N^2),fill(low_zero,N),fill(low_zero,N),ones(T,N),ones(T,N),fill(low_zero,N*mini_batch))
    upp = vcat(ones(T,N^2),p_bar,assets,Array{T}(fill(uppα,N)),Array{T}(fill(uppβ,N)),repeat(p_bar,mini_batch))

    ## linear constraints
    # hᵀz=q

    h₁₁ᵀ = hcat(sparse(kron(Matrix{T}(I,N,N),p_bar')))
    h₁ᵀ = hcat(h₁₁ᵀ,spzeros(T,N,N),sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
    q₁ = assets

    h₂₁ᵀ = repeat(spdiagm(p_bar),1,N)
    h₂ᵀ = hcat(h₂₁ᵀ,sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
    q₂ = p_bar

    hᵀ = vcat(h₁ᵀ,h₂ᵀ)
    q = vcat(q₁,q₂)
    function c_lin(z)
        return hᵀ*z .- q
    end

    ## quadratic constraints
    # zᵀHz = 0
    li = LinearIndices((N,N))
    liᵀ= transpose(li)
    id = li[tril!(trues(N,N), 0)]
    idᵀ = liᵀ[tril!(trues(N,N), 0)]
    H₃= sparse(vcat(id,idᵀ),vcat(idᵀ,id),T(0.5),M,M)
    h₃ = spzeros(T,M)
    q₃ = spzeros(T,1)

    H = H₃
    function c_quad(z)
        return dot(z,H*z) #transpose(z)*H*z
    end

    ## chance constraints ⟺ c_chance(z)=0
    # Prob(c[i]x[i]>=w[i])=delta[i] 


    beta_ccdf(a,b,x) = zero(x) < x < one(x) ? one(x) - IncBetaDer.beta_inc_grad(a, b, x)[1] : zero(x)
    function c_chance(z)
        zeroz = zero(z[1])
        onez = one(z[1])
        return beta_ccdf.(selα*z,selβ*z,max.(min.(w./(selc*z),onez),zeroz) ) .- delta
    end

    ## clearing vector constraint
    # quadratic version (ignores min{p_bar,max{0,...}} )
    function c_p(z,x)
        c = selc*z
        p = getp(z)
        Aᵀ = transpose(getA(z))
        return reshape((1+γ)*Aᵀ*p + ((1+γ)*(T(1) .- x)).*repeat(c,1,mini_batch)+γ*repeat(p_bar,1,mini_batch)-p,N*mini_batch)
    end

    # non-linear version
    function c_pNL(z,x)
        c = selc*z
        p = getp(z)
        Aᵀ = transpose(getA(z))
        return reshape(min.(max.(0,(1+γ)*Aᵀ*p + ((1+γ)*(T(1) .- x)).*repeat(c,1,mini_batch)+γ*repeat(p_bar,1,mini_batch)),p_bar)-p,N*mini_batch)
    end

    # write constraints as cL<=c(z)<=cU, setting cL=cU for equality constraints
    cL = spzeros(3N+1+mini_batch*N)
    cU = cL
    function c(z)
        #x = rand(N,mini_batch)
        x = ones(N,mini_batch)/10
        vcat(c_lin(z),
            c_quad(z),
            c_chance(z),
            c_p(z,x)
            )
    end

    # objective to max or min
    function obj(z)
        c = selc*z
        α = selα*z
        β = selβ*z
        p = getp(z)
        return sum(c.*α./(α + β))-sum(p)
    end

    # write constraints as cL<=c(z)<=cU, setting cL=cU for equality constraints
    # cL = spzeros(3N+1)
    # cU = cL
    # function c(z)
    #     vcat(hᵀ*z .- q,
    #         dot(z,H*z),
    #         beta_ccdf.(selα*z,selβ*z,w./(selc*z)) .- delta
    #         )
    # end
    # objective to max or min
    # function obj(z)
    #     Aᵀ = transpose(getA(z))
    #     c = selc*z
    #     α = selα*z
    #     β = selβ*z
    #     p = getp(z)
    #     d = Beta.(ones(N),2*ones(N))
    #     x = transpose(reduce(hcat,rand.(d, mini_batch)))
    #     return sum(c.*α./(α + β))-sum(min.(max.(0,(1+γ)*Aᵀ*p + ((1+γ)*(T(1) .- x)).*repeat(c,1,mini_batch)+γ*repeat(p_bar,1,mini_batch))))
    # end

    function obj_full(z)
        return obj(z)+sum(p_bar)
    end

    lowA = Array(low)
    uppA = Array(upp)

    z0 = (low+upp)/2

    m_chance = Nonconvex.Model(z->0.0)
    addvar!(m_chance,lowA,uppA)
    add_eq_constraint!(m_chance,z->c_chance(z))

    alg = IpoptAlg()
    options = Nonconvex.IpoptOptions()
    r_chance = Nonconvex.optimize(m_chance, alg, z0, options = options)

#     tol = 1e-16
#     @test r_chance.minimum ≈ 0 
#     @test sum(abs2,c_lin(r_chance.minimizer)) < tol

end