using ModelingToolkit
using Zygote
using ModelingToolkit.StructuralTransformations: numerical_nlsolve
using ChainRulesCore
using ForwardDiff

dfdu(f, u) = u isa Number ? ForwardDiff.derivative(f, u) : ForwardDiff.jacobian(f, u)
dfdp(f, u, p) = u isa Number && p isa Number ? ForwardDiff.derivative(f, p) :
                u isa Number ? ForwardDiff.gradient(f, p) :
                ForwardDiff.jacobian(f, p)

function ChainRulesCore.frule((_, _, _, ṗ), nlsolve::typeof(numerical_nlsolve), f, u0, p)
    u = nlsolve(f, u0, p)
    fu, fp = Base.Fix2(f, p), Base.Fix1(f, u)

    pushforward = -dfdu(fu, u) \ (dfdp(fp, u, p) * ṗ) # TODO: jvp
    return u, pushforward
end

function ChainRulesCore.rrule(nlsolve::typeof(numerical_nlsolve), f, u0, p)
    u = nlsolve(f, u0, p)
    fu, fp = Base.Fix2(f, p), Base.Fix1(f, u)

    function nlsolve_pullback(ū)
        p̄ = -dfdp(fp, u, p)' * (dfdu(fu, u)' \  ū) # TODO: vjp
        Zero(), Zero(), Zero(), p̄
    end
    return u, nlsolve_pullback
end

function nlfun(p)
    u0 = ones(2)
    numerical_nlsolve(u0, p) do u, p
        x, y = u
        [(x+3)*(y^3-7)+p[1]
         sin(y*exp(x)-p[2])]
    end |> sum
end

# Zygote.gradient(nlfun, [18, -1])[1]
ForwardDiff.gradient(nlfun, [18, -1])
#=
julia> Zygote.gradient(nlfun, [18, -1])[1]
2-element Vector{Float64}:
 -0.009532247816399934
  0.4895491994463628

julia> ForwardDiff.gradient(nlfun, [18, -1])
2-element Vector{Float64}:
 -0.009532385532365253
  0.48954879301835347
=#
