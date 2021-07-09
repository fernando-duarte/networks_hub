#include("NU2_preamble.jl")

# module NetworkModel

using NLPModels

# export NetNLP

mutable struct NetNLP <: AbstractNLPModel{Float64, Vector{Float64}}
    z0::AbstractVector
    low::AbstractVector
    upp::AbstractVector
    Q::Integer
    jac_sp::AbstractMatrix
    rows_jac_sp::AbstractVector
    cols_jac_sp::AbstractVector
    hess_sp::AbstractMatrix
    rows_hess_sp::AbstractVector
    cols_hess_sp::AbstractVector
    N::Integer
    minimize::Bool
    name::String
    p::AbstractVector                      #parameters
    meta::NLPModelMeta
    counters::Counters 
end

# show_header(io::IO, nlp::NetNLP) = println(io, "NetNLP - Network model")

# function Base.show(io::IO, nlp::NetNLP)
#   show_header(io, nlp)
#   show(io, nlp.meta)
#   show(io, nlp.counters)
# end

function NetNLP(z0,low,upp,Q,jac_sp,rows_jac_sp,cols_jac_sp,hess_sp,rows_hess_sp,cols_hess_sp,N; minimize = true,name = "Network Optimization",p=[])
  meta = NLPModelMeta(M; 
            x0 = z0,
            lvar = low,
            uvar = upp,
#             nlvb: number of nonlinear variables in both objectives and constraints
#             nlvo: number of nonlinear variables in objectives (includes nlvb)
#             nlvc: number of nonlinear variables in constraints (includes nlvb)
            nbv = 0,
            niv = 0,
            nlvbi = 0,
            nlvci = 0,
            nlvoi = 0,
            nwv = 0,
            ncon = Q,
#             y0: initial Lagrange multipliers
            lcon = zeros(Q),#zeros(ncon),
            ucon = zeros(Q),#zeros(ncon),
            #nnzo: number of nonzeros in all objectives gradients
            nnzj= nnz(jac_sp),
            nnzh = nnz(hess_sp),
            nlin = 2N,
            nnln = Q-2N,#ncon-nlin,
            nnnet = 0,
            nlnet = 0,
            lin = 1:2N,#1:nlin,
            nln = setdiff(1:Q,1:2N),#setdiff(1:ncon,1:nlin),
#             nnet: indices of nonlinear network constraints
#             lnet: indices of linear network constraints
            minimize = minimize,
            nlo = 1,
            islp = false,
            name = name
                )
  return NetNLP(z0,low,upp,Q,jac_sp,rows_jac_sp,cols_jac_sp,hess_sp,rows_hess_sp,cols_hess_sp, N,minimize,name,p,meta,Counters())
end

function NLPModels.obj(nlp::NetNLP, z::AbstractVector)
  return on(z,nlp.p)
end

function NLPModels.grad!(nlp::NetNLP, z::AbstractVector, g::AbstractVector)
  gn_iip(g,z,nlp.p)
  return g
end

function NLPModels.cons!(nlp::NetNLP, z::AbstractVector, c::AbstractVector)
  sn_iip(c,z,nlp.p)
  return c
end

function NLPModels.jac_coord!(
  nlp::NetNLP,
  z::AbstractVector,
  vals::AbstractVector
)
  jac_nzval_iip(vals,z)
  return vals
end

function NLPModels.jac_structure!(
  nlp::NetNLP,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  rows .= nlp.rows_jac_sp
  cols .= nlp.cols_jac_sp
  return rows, cols
end

function NLPModels.hess_structure!(
  nlp::NetNLP,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  rows .= nlp.rows_hess_sp
  cols .= nlp.cols_hess_sp
  return rows, cols
end

function NLPModels.hess_coord!(
  nlp::NetNLP,
  z::AbstractVector,
  vals::AbstractVector;
  obj_weight::Real = 1.0,
)
  hess_σ_nzval_iip(vals,z,obj_weight)
  return vals
end

function NLPModels.hess_coord!(
  nlp::NetNLP,
  z::AbstractVector,
  λ::AbstractVector,
  vals::AbstractVector;
  obj_weight::Real = 1.0,
)
  hess_λ_nzval_iip(vals,z,λ,obj_weight)
  return vals
end

function NLPModels.hprod!(
  nlp::NetNLP,
  z::AbstractVector,
  λ::AbstractVector,
  v::AbstractVector,
  Hv::AbstractVector;
  obj_weight::Real = 1.0,
)
  Hv_iip(Hv,z,λ,v,obj_weight)
  return Hv
end

# end # module