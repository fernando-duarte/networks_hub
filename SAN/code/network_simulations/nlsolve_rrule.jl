using NLsolve
using Zygote
using ChainRulesCore
using IterativeSolvers
using LinearMaps
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using Random
Random.seed!(1234)

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(nlsolve), f, x0; kwargs...)
    result = nlsolve(f, x0; kwargs...)
    function nlsolve_pullback(Δresult)
        Δx = Δresult[].zero
        x = result.zero
        _, f_pullback = rrule_via_ad(config, f, x)
        JT(v) = f_pullback(v)[2] # w.r.t. x
        # solve JT*Δfx = -Δx
        L = LinearMap(JT, length(x0))
        Δfx = gmres(L, -Δx)
        ∂f = f_pullback(Δfx)[1] # w.r.t. f itself (implicitly closed-over variables)
        return (NoTangent(), ∂f, ZeroTangent())
    end
    return result, nlsolve_pullback
end

const N = 10000
const nonlin = 0.1
const A = spdiagm(0 => fill(10.0, N), 1 => fill(-1.0, N-1), -1 => fill(-1.0, N-1))
const p0 = randn(N)
h(x, p) = A*x + nonlin*x.^2 - p
solve_x(p) = nlsolve(x -> h(x, p), zeros(N), method=:anderson, m=10).zero
obj(p) = sum(solve_x(p))

# need an rrule for h as Zygote otherwise densifies the sparse matrix A
# https://github.com/FluxML/Zygote.jl/issues/931
function ChainRulesCore.rrule(::typeof(h), x, p)
    y = h(x, p)
    function my_h_pullback(ȳ)
        ∂x = @thunk(A'ȳ + 2nonlin*x.*ȳ)
        ∂p = @thunk(-ȳ)
        return (NoTangent(), ∂x, ∂p)
    end
    return y, my_h_pullback
end

g_auto = Zygote.gradient(obj, p0)[1]
g_analytic = gmres((A + Diagonal(2*nonlin*solve_x(p0)))', ones(N))
display(g_auto)
display(g_analytic)
@show sum(abs, g_auto - g_analytic) / N  # 7.613631947123168e-17

@btime Zygote.gradient(obj, p0);                               # 17.825 ms (1288 allocations: 26.41 MiB)
@btime gmres((A + Diagonal(2*nonlin*solve_x(p0)))', ones(N));  # 18.214 ms (906 allocations: 24.03 MiB)
