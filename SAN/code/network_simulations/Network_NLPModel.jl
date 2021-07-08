include("NU2_preamble.jl")

module NetworkModel

export NetNLP

mutable struct NetNLP <: AbstractNLPModel
  p  #parameters
  meta::NLPModelMeta
  counters :: Counters 
end

show_header(io::IO, nlp::NetNLP) = println(io, "NetNLP - Network model")

function Base.show(io::IO, nlp::NetNLP)
  show_header(io, nlp)
  show(io, nlp.meta)
  show(io, nlp.inner.counters)
end

function NetNLP(p)
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
            lcon = zeros(ncon),
            ucon = zeros(ncon),
            #nnzo: number of nonzeros in all objectives gradients
            nnzj= nnz(jac_sp),
            nnzh = nnz(hess_sp),
            nlin = 2N,
            nnln = ncon-nlin,
            nnnet = 0,
            nlnet = 0,
            lin = 1:nlin,
            nln = nlin+1:nlin+ncon,
#             nnet: indices of nonlinear network constraints
#             lnet: indices of linear network constraints
            minimize = true,
            nlo = 1,
            islp = false,
            name = "Network Optimization"
                )
  return NetNLP(p, meta, Counters())
end

@default_counters NetNLP inner

function NLPModels.obj(nlp::NetNLP, z::AbstractVector)
  return on(z,nlp.p)
end

function NLPModels.grad!(nlp::NetNLP, z::AbstractVector, g::AbstractVector)
  gn_iip(g,z,nlp.p)
  return g
end

function NLPModels.hess_structure!(
  nlp::NetNLP,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  rows .= rows_hess_sp
  cols .= cols_hess_sp
  return rows, cols
end

function NLPModels.hess_coord!(
  nlp::NetNLP,
  z::AbstractVector,
  vals::AbstractVector;
  obj_weight::Real = 1.0,
)
  hess_σ_iip(vals,z,obj_weight)
  return vals
end

function NLPModels.hprod!(
  nlp::NetNLP,
  z::AbstractVector,
  v::AbstractVector,
  Hv::AbstractVector;
  obj_weight::Real = 1.0,
)
  hprod!(nlp.inner, x, v, Hv, obj_weight = obj_weight)
  Hv .+= nlp.ρ * obj_weight * v
  return Hv
end

end # module