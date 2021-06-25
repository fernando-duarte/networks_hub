import Pkg;Pkg.activate("nlp")
Pkg.instantiate()

using LinearAlgebra, DataFrames, SpecialFunctions,  SparseArrays, Random, BenchmarkTools, Distributions, JLD2, Missings, Test, NLsolve, MeasureTheory
using BSON, XLSX
using Parameters
using Dates, MonthlyDates
include("IncBetaDer.jl")

using Plots, Gadfly, Cairo, Fontconfig, LightGraphs, SimpleWeightedGraphs, GraphRecipes
using Setfield

# α=2.0; β=3.0;
# m= MeasureTheory.Beta(α, β)
# beta_logpdf(x) = logdensity(m, Lebesgue(𝕀), x) 

# define type Network
@with_kw struct Network{T,NN} @deftype Array{T}
    b = T[55.0,  55.0,  100.0,  55.0,  55.0] # outside liabilities
    c = T[50.0, 50.0, 150.0, 50.0, 50.0] # outside assets
    w = T[5.0, 5.0, 10.0, 5.0, 5.0]; # net worth
    p_bar = T[55.0, 55.0, 140.0, 55.0, 55.0]; # total liabilities
    delta = zeros(T,NN); # probabilities of default
    nm_b::Vector{Int64} = Int64[]
    nm_c::Vector{Int64} = Int64[]
    γ::T = T(0)
    dist::Product{Continuous} = Product(Distributions.Beta.(ones(NN),Inf)) # Distributions.jl only works with Float64 -- α==1.0 and β==Inf is consistent with delta.==0
    date::TimeType = today()
    name::Vector{String} = ["North-west node","Nort-east node","Central Node","South-east node","South-west node"]
    name_short::Vector{String} = ["NW","NE","CTR","SE","SW"]
    
    a = w .+ p_bar; # total assets
    d = a .- c; # inside assets
    f = p_bar .- b; # inside liabilities
    N::Integer = NN; # number of nodes
        
    A::Matrix{T} = T[0 0 0 0 0; 0 0 0 0 0; 10.0/p_bar[3] 10.0/p_bar[3] 0 10.0/p_bar[3] 10.0/p_bar[3]; 0 0 0 0 0; 0 0 0 0 0] # matrix of relative liabilities
end

## constructors   
function netEmp(df::DataFrame,
            N::Integer,
            A::Union{Nothing,Matrix}=nothing,
            dist::Union{Nothing,AbstractMeasure,Sampleable} = nothing;
            T::Type=Float64,
            dates::Vector{Date}=[Date(QuarterlyDate("2008-Q4"))],
            γ::Real=0, 
        )
    df.julia_dates = toJuliaDates(df.qt_dt)
    net = Array{Network{T,N},1}(undef, length(dates))
    for i=1:length(dates)
        keep = df.julia_dates.==lastdayofquarter(dates[i])
        df_t = sort(df[keep,:], :assets, rev = true)
  
        nm_b = [findall(x->x==false,ismissing.(df_t.b))...]
        nm_b = nm_b[isless.(nm_b, N+1)] # indices for non-missing b
        nm_c = [findall(x->x==false,ismissing.(df_t.c))...]
        nm_c = nm_c[isless.(nm_c, N+1)] # indices for non-missing b
        df_t = coalesce.(df_t, df_t.assets/1.5) # replace missing with some value
        
        b = [deepcopy(Array{T}(df_t[1:N,[:b]]))...]
        c = [deepcopy(Array{T}(df_t[1:N,[:c]]))...]
        w = [deepcopy(Array{T}(df_t[1:N,[:w]]))...]
        p_bar = [deepcopy(Array{T}(df_t[1:N,[:p_bar]]))...]
        delta = [deepcopy(Array{T}(df_t[1:N,[:delta]]))...] 
               
        # create A consistent with rest of network
        if isnothing(A)
            a = w .+ p_bar; # total assets  
            mat_ind = randIndA(N) # random edges
            A = zeros(T,N,N)
            A[mat_ind] .= solve_constrA(mat_ind,b,p_bar,c,a)
        end

        # make α,β parameters of the beta distribution consistent with prob of default
        if isnothing(dist)
            α = ones(N)
            β = solve_constrβ(α,w,c,delta)    
            β[delta.==0] .= Inf
            dist = Product(Distributions.Beta.(α,β)) 
        end
        net[i]= Network{T,N}(      
            b = b,
            c = c,
            w = w,
            p_bar = p_bar,
            delta = delta, 
            nm_b = nm_b,
            nm_c = nm_c,
            γ=γ,
            date=dates[i],
            name = [deepcopy(Array(df_t[1:N,[:nm_short]]))...], 
            name_short = [deepcopy(Array(df_t[1:N,[:tkr]]))...], 
            A = A,
            dist = dist,
        )
    end #for
    return net
end
 
function netEmp(fileName::String,N::Integer;
            T::Type=Float64,
            dates::Vector{Date}=[Date(QuarterlyDate("2008-Q4"))],
            γ::Real=0,
            dist::Product{Continuous} = Product(Distributions.Beta.(ones(T,N),T[2.0]))
            )         
        #create data::DataFrame from excel file called fileName
        xf = XLSX.readxlsx(fileName) 
        data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) 
        unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
        # convert stata %tq date format into Julia type QuarterlyDate <: TimeType
        data.julia_dates = toJuliaDates(convert.(Int,data.qt_dt))
        scale = floor(Int,log10(maximum(skipmissing(data[:,:assets]))))
        units = 10^(scale);
        data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units
        # create network(s)
        data[:,:delta] .= any(names(data).=="delta") ? data[!,:delta] : zeros(T,size(data,1))
        netEmp(data,N;T=T,dates=dates,γ=γ)          
end

## plot recipe
function plot(network::Network;plot_outside::Bool=false)
    network_round = roundNet(network)
    @unpack b,c,w,p_bar,N,A = network_round
    A[A.<1e-3] .= 0   # ignore small edges
    #A .= round.(A,digits=2)
    
    
    if plot_outside    
        # create a relative liabilities matrix that includes extra nodes for outside assets and liabilities. For each node inside the network, we define two extra nodes, one for outside assets and another for outside liabilities.      
        Aall = vcat(hcat(p_bar.*A,Diagonal([c...]),zeros(N,N)),hcat(vcat(zeros(N,N),Diagonal([b...])),zeros(2N,2N))) # relative liabilities matrix
        Aall[Aall.<1e-3] .= 0
        Nall = size(Aall,1) # number of nodes

        p=graphplot(LightGraphs.DiGraph(Aall),
                  nodeshape= :circle,
                  markersize = 0.2,
                  node_weights =  [hcat([w...],[b...],[c...]).^(1/30)...],
                  markercolor = :white, #range(colorant"yellow", stop=colorant"red", length=N+1),
                  names =  hcat(w...,b...,c...),
                  fontsize = 8,
                  edgecolor = :black, # for GY example: hcat(:red,:red,:black,:black,:black,:black,:red,:red,:red,:red,:red,:red,:red,:red),
                  linestyle = :solid, # for GY example: hcat(:dashdot,:dashdot,:solid,:solid,:solid,:solid,:dashdot,:dashdot,:dashdot,:dashdot,:dashdot,:dashdot,:dashdot,:dashdot),
                  edgewidth = (s,d,w)->0.02*Aall[s,d],
                  arrow = true,
                  method = :stress,
                  curves = false,
                  edgelabel = Aall,
                  axis_buffer = 0.25,
                  edgelabel_offset = 0,
                  shorten = 0,
                  edge_label_box = true,
                  layout_kw = Dict(:iterations => 10000),
                  size = (1000,1000),
        )
    else
        p=graphplot(LightGraphs.DiGraph(A),
              nodeshape = :circle,
              markersize = 0.1,
              node_weights = w,
              markercolor = :white, #range(colorant"yellow", stop=colorant"red", length=N+1),
              names = w,
              fontsize = 8,
              linecolor = :darkgrey,
              edgewidth = (s,d,w)->15*A[s,d], # width proportional to relative liabilities
              arrow = true,
              method = :sfdp,
              curves = false,
              edgelabel = round.(p_bar.*A,digits=2),
        ) 
    # attributes here: https://docs.juliaplots.org/latest/generated/graph_attributes/
    # method `:spectral`, `:sfdp`, `:circular`, `:shell`, `:stress`, `:spring`, `:tree`, `:buchheim`, `:arcdiagram` or `:chorddiagram`.
    end
end    
    
## utility
function toJuliaDates(stata_tq::Vector{T}) where T
    if T==Date
        return stata_tq        
    else T==Int
        # convert stata %tq date format into Julia type QuarterlyDate <: TimeType
        EPOCH = Dates.Date(1960, 1, 1)
        elapsed_quarters(Δ::Integer) = EPOCH + Dates.Quarter(Δ)
        stata_quarters_to_date(q::Integer) = lastdayofquarter(elapsed_quarters(q))        
        return stata_quarters_to_date.(stata_tq)
    end
    return nothing        
end

function randIndA(N)    
    lt = LowerTriangular(ones(N,N))-Matrix(I,N,N) 
    ind = findall( lt .> 0 )
    indᵀ = map(x->CartesianIndex(x[2],x[1]),CartesianIndices(lt)[ind])
    d = Distributions.Binomial(1,1/2)
    v = rand(d,length(ind))
    return [v[i]*LinearIndices(lt)[ind[i]] + (1-v[i])*LinearIndices(lt)[indᵀ[i]] for i=1: length(ind)]
end
      
function constrA(x::Array{T},mat_ind::Vector{Int},b::Vector{R},p_bar::Vector{R},c::Vector{R},a::Vector{R}) where {T<:Real,R<:Real}
    N = length(p_bar);
    A = zeros(T,N,N); 
    A[mat_ind] .= x
    return vcat(b - p_bar + diagm(p_bar)*A*ones(N), c - a + transpose(A)*p_bar)
end

function solve_constrA(mat_ind::Vector{Int},b::Vector{R},p_bar::Vector{R},c::Vector{R},a::Vector{R}) where {R<:Real}
  constrA_x = (x) -> constrA(x,mat_ind,b,p_bar,c,a) # this creates a new function of one variable (x) which "closes over" the current values of `mat_ind,b,p_bar,c,a`
  x0 = ones(length(mat_ind))
  nlsolve(constrA_x,x0,autodiff = :forward).zero
end

function constrβ(β::Array{T},α::Vector{R},w::Vector{R},c::Vector{R},delta::Vector{R}) where {T<:Real,R<:Real}
    # cdf(Beta) only accepts Float64
   #1 .- cdf.(Beta.(α,β),w ./ c) .- delta
   beta_ccdf(a,b,x) = zero(x) < x < one(x) ? one(x) - IncBetaDer.beta_inc_grad(a, b, x)[1] : zero(x)
   beta_ccdf.(α,β,w ./ c) .- delta     
end

function solve_constrβ(α::Vector{R},w::Vector{R},c::Vector{R},delta::Vector{R}) where {R<:Real}
  constrβ_x = (x) -> constrβ(x,α,w,c,delta)
  β0 = fill(30.0,length(α))
  nlsolve(constrβ_x,β0,autodiff = :forward).zero
end

function roundNet(net::Network)
    # rounds numbers inside properties of net
    pn = propertynames(net)
    reals_ind = ([eltype(getproperty(net,p)) for p in pn] .<: Real)
    reals_name = pn[reals_ind]
    reals_name = reals_name[findall(x->x!=:dist,reals_name)]
    reals_val = [round.(getproperty(net,p),digits=2) for p in reals_name if p!=:dist]
    notreals_names = vcat(setdiff(pn,reals_name),:dist)
    noreals_val = [getproperty(net,p) for p in notreals_names]
    nm = vcat(reals_name...,notreals_names)
    val = vcat(reals_val,noreals_val)
    function copy_with_modification(original::T, field_to_change, new_value) where {T}
        val(field) = field==field_to_change ? new_value : getfield(original, field)
        T(val.(fieldnames(T))...)
    end
    for p in zip(nm,val)
        net = copy_with_modification(net,p[1],p[2])
    end
    return net
end


## instantiate
    
# empty constructor defaults to an example from Glasserman and Young
T=Float64
N=5
netGY = Network{T,5}()
    
# network from xlsx file
fileName="node_stats_forsimulation_all.xlsx"
N=5
net1 = netEmp(fileName,N)[1]
    
fileName="glasserman_young_example.xlsx"
N = 5 
netGY2 = netEmp(fileName,N)[1]
 

# network from DataFrame
data_dict = BSON.load("data.bson")
net2 = netEmp(data_dict[:data],5)[1]

# plots
pGY = plot(netGY)
pGY_in = plot(netGY,plot_outside = false)
pGY_out = plot(netGY,plot_outside = true)
savefig(pGY, "GYnet.pdf") 
savefig(pGY_in, "GYnet_in.pdf") 
savefig(pGY_out, "GYnet_out.pdf") 

pGY2 = plot(netGY2)
savefig(pGY2, "GYnet.pdf") 

pnet1=plot(net1)
savefig(pnet1, "net.pdf") 

## optimization
       
#= Notes
f(betas) = h(betas...,1,1,1,1,1)

Optim.optimize(
    f,
    ones(6), 
    method = ConjugateGradient(),
    autodiff = :forward,
    x_tol = 1e-7
)
        
        
=#        
mini_batch = 1;
M = N^2+4*N+N*mini_batch

# z = vcat(A...,b,c,α,β,p...)
# p is N-by-mini_batch

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
lowα = 0.5 ;uppα = 10
lowβ = 1   ;uppβ = 150
low_zero = eps(T) # eps(T) T(0.005)

low = vcat(spzeros(T,N^2),fill(low_zero,N),fill(low_zero,N),fill(lowα,N),fill(lowβ,N),fill(low_zero,N*mini_batch))
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
    zeroz = zero(z[1])+0.0001
    onez = one(z[1])-0.0001
#     display(selc*z)
    return beta_ccdf.(selα*z,selβ*z,max.(min.(w./(selc*z),onez),zeroz) ) .- delta
end

## clearing vector constraint
# quadratic version (ignores min{p_bar,max{0,...}} )
function c_p(z,x)
    c = selc*z
    p = getp(z)
    Aᵀ = transpose(getA(z))
    return reshape(p - ((1+γ)*Aᵀ*p + ((1+γ)*(T(1) .- x)).*repeat(c,1,size(x,2))-γ*repeat(p_bar,1,size(x,2))) ,N*size(x,2))
end

# non-linear version
function c_pNL(z,x)
    c = selc*z
    p = getp(z)
    Aᵀ = transpose(getA(z))
    p_bar_rep = repeat(p_bar,1,size(x,2))
    return reshape(p-min.(max.(0,(1+γ)*Aᵀ*p + ((1+γ)*(T(1) .- x)).*repeat(c,1,size(x,2))-γ*p_bar_rep),p_bar_rep),N*size(x,2))
end

# write constraints as cL<=c(z)<=cU, setting cL=cU for equality constraints
cL = spzeros(3N+1+mini_batch*N)
cU = cL
function c_(z,x)
    #x = rand(N,mini_batch)
    #x = ones(N,mini_batch)/10
    vcat(c_lin(z),
        c_quad(z),
        c_chance(z),
        c_p(z,x)
        )
end

function c_NL(z,x)
    #x = rand(N,mini_batch)
    #x = ones(N,mini_batch)/10
    vcat(c_lin(z),
        c_quad(z),
        c_chance(z),
        c_pNL(z,x)
        )
end

function c_(z)
    vcat(c_lin(z),
        c_quad(z),
        c_chance(z)
        )
end

function c_NL(z)
    vcat(c_lin(z),
        c_quad(z),
        c_chance(z)
        )
end

# objective to max or min
function obj(z)
    c = selc*z
    α = selα*z
    β = selβ*z
    p = getp(z)
    return sum(c.*α./(α + β))-sum(mean(p,dims=2))
end

function obj(z,x)
    Aᵀ = transpose(getA(z))
    c = selc*z
    α = selα*z
    β = selβ*z
    p_guess = getp(z) 
    p_guess_rep = repeat(p_guess,1,size(x,2))
    p_bar_rep = repeat(p_bar,1,size(x,2))
    function clearing_p!(F, p, x, Aᵀ, c)
        crep = repeat(c,1,size(x,2))
        F .= p.-min.(p_bar_rep,max.( (1+γ)*(Aᵀ*p .+ (one(x[1]) .- x).*crep) .- γ.*p_bar_rep, 0 ) )
    end
    p = nlsolve((F, p)->clearing_p!(F, p, x, Aᵀ, c),p_guess_rep, autodiff = :forward)
    #@test converged(p)
    #return p.zero
    return sum(c.*x)-sum(mean(p.zero,dims=2)) # spillovers for network 
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end

function obj_full(z)
    obj(z)+sum(p_bar)
end

function obj_full(z,x)
    obj(z,x)+sum(p_bar)
end

function loss_d(z,x)
    cx = c.*x
    return sum(xc+max.(0,cx-w)) # spillovers for disconnected network 
end

end # module