module NetworkType

using Convex, SCS, NLsolve
using LinearAlgebra, Distributions, Random, SpecialFunctions,  SparseArrays, MeasureTheory
using BenchmarkTools, Test
using DataFrames, BSON, XLSX, JLD2, Missings
using Dates, MonthlyDates
using Plots, Gadfly, Cairo, Fontconfig, LightGraphs, SimpleWeightedGraphs, GraphRecipes
using Setfield, Parameters, Suppressor
using MathOptInterface
const MOI = MathOptInterface
include("IncBetaDer.jl")

export Network, netEmp, plot, constrA, solve_constrA, constrβ, solve_constrβ, roundNet

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




##############################################################################################
##############################################################################################
##############################################################################################
export getA, trainingSetA, randIndA, allIndA

"""
    getA(net [, nonzero_ind, con])

Draw a random relative liabilities matrix A for a network `net` that has 
elements between 0 and 1 and satisfies
the constraints specified by `con` with possible values:

- `"diag"`: diagonal is zeros(N)
- `"netting"`: A[i,j]*A[j,i]==0 for i,j=1:N
- `"liabs"`: b - p_bar + diagm(p_bar)*A*ones(N) == 0, ie, payments to other nodes add up to inside liabilities 
- `"assets"`: c - a + transpose(A)*p_bar == 0, ie payments from other nodes add up to inside assets 
- `"all"`: satisfies "netting", "liabs" and "assets"

To satisfy a combination of constraints, `con` can be an array, e.g. `con=["netting","liabs"]`. 

# Arguments
- `net`: Either a struct with fields `b`, `c`, `p_bar` and `a`, or an integer `N`. If `net` is an integer, 
          `b`, `c`, `p_bar` and `a` must be in the scope of `getA` (e.g., be globals)
# Optional kwargs
- `con::Union{String,Vector{String}}="all"`: constraints that A satisfies
- `nonzero_ind::Union{Nothing,Vector{Int}}=nothing`: linear index for elements of A allowed to be non-zero

"""
function getA(net;nonzero_ind::Union{Nothing,Vector{Int},Vector{Vector{Int64}}}=nothing,con::String="all",optimizer=SCS.Optimizer(verbose=false, eps=1e-12,linear_solver=SCS.DirectSolver))
    # if nonzero_ind is a vector of vectors, process one a time
    if typeof(nonzero_ind)==Vector{Vector{Int64}}
        Q = length(nonzero_ind)
        sampleA = Array{Union{Any,Nothing},1}(undef,Q) # initialize
        for i=1:Q
            sampleA[i] = getA(net;nonzero_ind=nonzero_ind[i],con=con)   
        end
        return sampleA
    end
      
    # if nonzero_ind is a vector
    if isa(net,Integer) 
        N=net 
    else 
        @unpack b,c,p_bar,a,N,γ = net 
    end
    if isa(con,String) con = [con] end
    if con==["all"] con = ["netting","liabs","assets"] end
    fc(s) = any(map(x->x==s, con))
    if (fc("netting") && fc("diag")) filter!(x->x!="diag",con) end
    
    if con == ["diag"] 
        A = rand(N,N)
        A = A - diagm(diag(A))
        return A
    end
    if (isnothing(nonzero_ind) && fc("netting")) nonzero_ind = randIndA(N) end
    if con == ["netting"] 
        A = zeros(N,N)
        A[nonzero_ind] = rand(length(nonzero_ind))
        return A
    end
     
    #x = rand(N)
    zeros_ind = setdiff(1:N^2,nonzero_ind)
    A = zeros(N,N)
    
    AA = Variable(N,N, Positive()) 
    bb = Variable(N,Positive()) 
    cc = Variable(N,Positive()) 
    
    #pp = Variable(N,Positive())
    #fix!(bb,b) # free!(bb)
    #fix!(cc,c) # free!(cc)
    #fix!(pp,p_bar) # free!(pp)
    #problem = minimize(sum(pp))
    problem = minimize(0.0)
    problem.constraints = [ AA >= 0.0, AA <= 1.0 ,
                            bb >= 0.0, bb <= p_bar ,
                            cc >= 0.0, cc <= a,
                            #pp >= 0.0, pp <= p_bar                    
                            ];
    if fc("diag") push!(problem.constraints, diag(AA) == 0) end # A[i,i]==0
    if fc("netting") push!(problem.constraints, AA[zeros_ind] == 0) end # A[i,j]*A[j,i]==0
    if fc("liabs") push!(problem.constraints, bb - p_bar + diagm(p_bar)*AA*ones(N) == 0)  end # payments to other nodes add up to inside liabilities 
    if fc("assets") push!(problem.constraints, cc - a + transpose(AA)*p_bar == 0) end # payments from other nodes add up to inside assets 
    #if fc("p(x)") push!(problem.constraints, pp == (1+γ)*(transpose(AA)*pp + cc - diagm(x)*cc) - γ*p_bar) end # p(x) is a clearing vector
    

    @suppress_err begin
        Convex.solve!(problem,  optimizer, warmstart=false; silent_solver=true)
    end
    if problem.status == MOI.OPTIMAL #!(problem.status == MOI.INFEASIBLE || problem.status == MOI.ALMOST_INFEASIBLE || problem.status == MOI.INFEASIBLE_OR_UNBOUNDED)
        return (evaluate(AA),evaluate(bb),evaluate(cc))
    else
        return nothing
    end
end

# if nonzero_ind is a vector of vectors, process one a time
# function getAllA(net;nonzero_ind::Vector{Vector{Int64}},con::String="all")
#     Q = length(nonzero_ind)
#     sampleA = Array{Array{Float64,2},1}(undef,Q) # initialize
#     for i=1:Q
#         sampleA[i] = getA(net;nonzero_ind=nonzero_ind[i],con=con)   
#     end
#     return sampleA
# end

# randomly generate indices for A where it is allowed to be non-zero
function randIndA(N::Integer)    
    lt = LowerTriangular(ones(N,N))-Matrix(I,N,N) 
    ind = findall( lt .> 0 )
    indᵀ = map(x->CartesianIndex(x[2],x[1]),CartesianIndices(lt)[ind])
    d = Distributions.Binomial(1,1/2)
    v = rand(d,length(ind))
    return [v[i]*LinearIndices(lt)[ind[i]] + (1-v[i])*LinearIndices(lt)[indᵀ[i]] for i=1: length(ind)]
end

# generate all possible combinations of indices such that A[i,j]*A[j,i]==0
function allIndA(N::Integer)    
    lt = LowerTriangular(ones(N,N))-Matrix(I,N,N) 
    ind = findall( lt .> 0 )
    indᵀ = map(x->CartesianIndex(x[2],x[1]),CartesianIndices(lt)[ind])
    M = length(ind)
    ind01 = collect.(reverse.(Iterators.product(fill([true,false],M)...))[:])
    #all_ind01 = [ind[ind01[i]] for i=1:length(ind01)]
    all_ind01 = [vcat(ind[ind01[i]],indᵀ[.!ind01[i]]) for i=1:length(ind01)]
    all_ind01_lin = [LinearIndices(lt)[all_ind01[i]] for i=1:length(all_ind01)]
    lin_ind = LinearIndices(lt)[ind];
    return (all_ind01_lin,all_ind01,ind01,lin_ind,ind)
end


    
function trainingSetA(net; sample_size::Integer = 20, max_iter::Integer = 2000, nonzero_ind::Union{Nothing,Vector{Int}}=nothing,con::String="all")
    if isa(net,Integer) 
        N=net 
    else 
        @unpack b,c,p_bar,a,N = net 
    end
    if max_iter<sample_size
       max_iter = sample_size 
    end
    q = 0 # number of matrices found after max_iter iterations
    sampleA = Array{Array{Float64,2},1}(undef,sample_size) # initialize
    iter = 1 
    while iter<=max_iter && q<sample_size
        newA = getA(net, nonzero_ind=nonzero_ind, con=con)
        if !isnothing(newA)
            sampleA[q+1] = newA
            q+=1
        end
        iter+=1
    end
                    
    if q==0
        @warn "Problem may be infeasible"
        return nothing                
    else
        return sampleA[1:q]
    end
end


end #module