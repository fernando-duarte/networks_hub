module generateA

using Convex, SCS, LinearAlgebra, Distributions, Test, Random, MathOptInterface, Parameters, Suppressor
const MOI = MathOptInterface

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
function getA(net;nonzero_ind::Union{Nothing,Vector{Int},Vector{Vector{Int64}}}=nothing,con::String="all")
    # if nonzero_ind is a vector of vectors, process one a time
    if typeof(nonzero_ind)==Vector{Vector{Int64}}
        Q = length(nonzero_ind)
        sampleA = Array{Union{Array{Float64,2},Nothing},1}(undef,Q) # initialize
        for i=1:Q
            sampleA[i] = getA(net;nonzero_ind=nonzero_ind[i],con=con)   
        end
        return sampleA
    end
      
    # if nonzero_ind is a vector
    if isa(net,Integer) 
        N=net 
    else 
        @unpack b,c,p_bar,a,N = net 
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
     
    zeros_ind = setdiff(1:N^2,nonzero_ind)
    A = Variable(N,N,Positive()) 
    problem = minimize(0.0)
    problem.constraints = [ A >= 0.0, A <= 1.0 ];
    if fc("diag") push!(problem.constraints, diag(A) == 0) end # A[i,i]==0
    if fc("netting") push!(problem.constraints, A[zeros_ind] == 0) end # A[i,j]*A[j,i]==0
    if fc("liabs") push!(problem.constraints, b - p_bar + diagm(p_bar)*A*ones(N) == 0)  end # payments to other nodes add up to inside liabilities 
    if fc("assets") push!(problem.constraints, c - a + transpose(A)*p_bar == 0) end # payments from other nodes add up to inside assets 
    @suppress_err begin
        Convex.solve!(problem,  SCS.Optimizer(verbose=false, eps=1e-12,linear_solver=SCS.DirectSolver), warmstart=false; silent_solver=true)
    end
    if !(problem.status == MOI.INFEASIBLE || problem.status == MOI.ALMOST_INFEASIBLE)
        return evaluate(A)
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

end # module


