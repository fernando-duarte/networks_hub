# I should really extend the methods in SpecialFunctions, but I have no clue
# why they return a tuple...
beta_inc_(p::Float64, q::Float64, x::Float64) = beta_inc(p, q, x)[1]
beta_inc_(p, q, x) = beta_inc_(promote(p, q, x)...)

# Since we will be having one specific use case (for now), we may be happy with
# implementing the specific dispatch we need...
# quantile(Beta(a,b), x) == beta_inc_inv(a, b, x)[1]
function beta_inc_(p::Dual{T}, q::Dual{T}, x::Dual{T}) where T
    z, dp, dq, dx = beta_inc_grad(value(p), value(q), value(x))
    d = [dp, dq, dx]
    input_parts = [partials(p), partials(q), partials(x)]
    input_parts_ = zip(collect.(input_parts)...)
    #┌ Info: these are the input partials wrapped up
    #│   collect(input_parts_) =
    #│    3-element Array{Tuple{Float64,Float64,Float64},1}:
    #│     (1.0, 0.0, 0.129359118670396)
    #│     (0.0, 1.0, -0.055141152997250496)
    #└     (0.0, 0.0, 0.6326024899526733)
    parts = map(x->sum(d .* x), input_parts_) 
    parts = ForwardDiff.Partials((parts...,))
    Dual{T}(z, parts) 
end

# The method is based on the identity: 
#   I(x;p,q) = K(x;p,q) * ₂F₁(1-q, 1; p+1; -x/(1-x)) 
# here, K(x;p,q) is the following function:
Kfun(x, p, q) = (x^p * (1. - x)^(q-1)) / (p * beta(p, q))

# The method uses a continued fraction representation to evaluate the
# hypergeometric series in the identity above.  The following f appears in the
# CF representation.
ffun(x, p, q) = q*x/(p*(1. - x))

# The {aₙ} sequence of the continued fraction for ₂F₁(1-q, 1; p+1; -x/(1-x)) 
a1fun(p, q, f) = p*f*(q - 1.)/(q*(p + 1.))
anfun(p, q, f, n) = n == 1 ? a1fun(p, q, f) : 
    p^2 * f^2 * (n - 1.) * (p + q + n - 2.) * (p + n - 1.) * (q - n) / 
        (q^2 * (p + 2n - 3.) * (p + 2n - 2.)^2 * (p + 2n - 1.))

# ... and the {bₙ} sequence.
function bnfun(p, q, f, n) 
    x = 2*(p*f + 2q)*n^2 + 2*(p*f + 2q)*(p - 1.)*n + p*q*(p - 2. - p*f) 
    y = (q * (p + 2n - 2.) * (p + 2n))
    x/y
end

# The method relies on the differentiating through the continued fraction. 
# We will need the derivatives of K(x;p, q) and the {aₙ} and {bₙ}
# The following is all straight out of the appendix of Boik & Robinson-Cox
dK_dp(x, p, q, K, ψpq, ψp) = K*(log(x) - 1.0/p + ψpq - ψp) 
dK_dq(x, p, q, K, ψpq, ψq) = K*(log(1.0 - x) + ψpq - ψq) 
function dK_dpdq(x, p, q)
    ψ = digamma(p+q)
    Kf = Kfun(x, p, q)
    dKdp = dK_dp(x, p, q, Kf, ψ, digamma(p))
    dKdq = dK_dq(x, p, q, Kf, ψ, digamma(q))
    dKdp, dKdq
end

da1_dp(p, q, f) = -p * f * (q - 1.) / (q * (p + 1.)^2)
function dan_dp(p, q, f, n)
    n == 1 && return da1_dp(p, q, f) 
    x = -(n - 1.) * f^2 * p^2 * (q - n)
    y = (-1. + p + q)*8n^3 + 
        (16p^2 + (-44. + 20q)*p + 26. - 24q)*n^2 + 
        (10p^3 + (14q - 46.)*p^2 + (-40q + 66.)*p - 28. + 24q)*n + 
         2p^4 + (-13 + 3q)*p^3 + (-14 + 30)*p^2 + 
        (-29 + 19q)*p + 10. - 8q
    z = q^2 * (p + 2n - 3.)^2 * (p + 2n - 2.)^3 * (p + 2n - 1.)^2
    x * y / z
end

da1_dq(p, q, f) = f * p / (q * (p + 1.))
function dan_dq(p, q, f, n)
    n == 1 && return da1_dq(p, q, f)
    x = p^2 * f^2 * (n - 1.) * (p + n -1.) * (2q + p - 2.)
    y = q^2 * (p + 2n - 3) * (p + 2n - 2)^2 * (p + 2n - 1.)
    x/y
end

function dbn_dp(p, q, f, n)
    x = (1. - p - q) * 4n^2 + (4p - 4. + 4q - 2p^2)*n + p^2 * q
    y = q * (p + 2n - 2.)^2 * (p + 2n)^2
    p * f * x / y 
end

dbn_dq(p, q, f, n) = - p^2 * f / (q * (p + 2n - 2.) * (p + 2n))

"""
    beta_inc_grad(a, b, x, [maxapp, minapp, ϵ])

Compute the partial derivatives of the incomplete Beta function with respect
to `a` and `b`. Returns the following tuple:
```math
(I_{x}(a, b), \\frac{\\partial I_{x}(a, b)}{\\partial a},
    \\frac{\\partial I_{x}(a, b)}{\\partial b})
```
This uses the method described in Boik & Robinson-Cox (1998), which uses
a continued fraction representation of the incomplete Beta function and its
partial derivatives.

The `maxapp` and `minapp` arguments specify the maximum resp. minimum number of
approximants to use in the continued fraction evaluation, while `ϵ` determines
the convergence threshold.
"""
function beta_inc_grad(a, b, x, maxapp=200, minapp=3, ϵ=1e-12)
    # Here I use `a` and `b` in the args for consistency with `beta_inc` 
    # However in the function body a ≡ p and b ≡ q.
    # x == 1 or 0 cases
    if x == one(x)
        return one(x), zero(x), zero(x), zero(x)
    elseif x == zero(x)
        return zero(x), zero(x), zero(x), zero(x)
    end
    dx = x^(a - 1.) * (1. - x)^(b - 1.) / beta(a,b)  # partial derivative with respect to x
    # swap tails if necessary
    if x <= a / (a + b)
        p = a
        q = b
        swap = false
    else
        x = 1. - x
        p = b
        q = a
        swap = true
    end
    # compute once
    K = Kfun(x, p, q)
    dK_dp, dK_dq = dK_dpdq(x, p, q)
    f = ffun(x, p, q)
    App = 1.
    Ap  = 1.
    Bpp = 0. 
    Bp  = 1. 
    dApp_dp = 0.
    dBpp_dp = 0.
    dAp_dp  = 0.
    dBp_dp  = 0.
    dApp_dq = 0.
    dBpp_dq = 0.
    dAp_dq  = 0.
    dBp_dq  = 0.
    dI_dp   = NaN
    dI_dq   = NaN
    Ixpq    = NaN
    Ixpqn   = NaN
    for n=1:maxapp
        An, Bn, an, bn = _nextapp(f, p, q, n, App, Ap, Bpp, Bp)
        
        # update derivatives for p
        dan = dan_dp(p, q, f, n)
        dbn = dbn_dp(p, q, f, n)
        dAn_dp = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dp, dAp_dp)
        dBn_dp = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dp, dBp_dp)
        
        # update derivatives for q
        dan = dan_dq(p, q, f, n)
        dbn = dbn_dq(p, q, f, n)
        dAn_dq = _dnextapp(an, bn, dan, dbn, App, Ap, dApp_dq, dAp_dq)
        dBn_dq = _dnextapp(an, bn, dan, dbn, Bpp, Bp, dBpp_dq, dBp_dq)

        # compute target
        Cn = An/Bn
        dI_dp = dK_dp * Cn + K * ((1. / Bn) * dAn_dp - (An/Bn^2) * dBn_dp)
        dI_dq = dK_dq * Cn + K * ((1. / Bn) * dAn_dq - (An/Bn^2) * dBn_dq)
        Ixpqn = K * Cn 

        # check convergence
        #@info "" K f dan dbn dAn_dq dBn_dq An Bn
        (abs(Ixpqn - Ixpq) < ϵ && n >= minapp) && break

        # update
        Ixpq = Ixpqn 
        App  = Ap
        Bpp  = Bp 
        Ap   = An
        Bp   = Bn
        dApp_dp = dAp_dp
        dApp_dq = dAp_dq
        dBpp_dp = dBp_dp
        dBpp_dq = dBp_dq
        dAp_dp  = dAn_dp
        dAp_dq  = dAn_dq
        dBp_dp  = dBn_dp
        dBp_dq  = dBn_dq
    end         
    if swap
        return 1. - Ixpqn, -dI_dq, -dI_dp, dx
    else
        return Ixpqn, dI_dp, dI_dq, dx
    end
end

# get the next approximant
function _nextapp(f, p, q, n, App, Ap, Bpp, Bp) 
    an = anfun(p, q, f, n)
    bn = bnfun(p, q, f, n)
    An = an*App + bn*Ap
    Bn = an*Bpp + bn*Bp
    An, Bn, an, bn
end

function _dnextapp(an, bn, dan, dbn, Xpp, Xp, dXpp, dXp)
    dan * Xpp + an * dXpp + dbn * Xp + bn * dXp 
end

# For testing, recursive evaluation of continued fraction.
function recursive_cf(x, p, q, depth)
    K = Kfun(x, p, q)
    f = ffun(x, p, q)
    function cf(n)
        n == depth && return bnfun(p, q, f, n)
        anfun(p, q, f, n)/(bnfun(p, q, f, n) + cf(n+1))
    end
    K*(1. + cf(1))
end

