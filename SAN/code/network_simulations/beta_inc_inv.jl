# These are copied from SpecialFunctions, but allowing other Real types to pass
# through the numerical methods (for autodiff)

beta_inc_inv_(a, b, p; lb = logbeta(a,b)[1]) =
    beta_inc_inv_(promote(a, b, p)..., lb=lb)

function beta_inc_inv_(a::T, b::T, p::T; lb=logbeta(a,b)[1]) where T
    q = 1. - p
    fpu = 1e-30
    x = p
    if p == 0.0
        return 0.0
    elseif p == 1.0
        return 1.0
    end

    #change tail if necessary

    if p > 0.5
        pp = q
        aa = b
        bb = a
        indx = true
    else
        pp = p
        aa = a
        bb = b
        indx = false
    end

    #Initial approx

    r = sqrt(-log(pp^2))
    pp_approx = r - @horner(r, 2.30753e+00, 0.27061e+00) / @horner(r, 1.0, .99229e+00, .04481e+00)

    if a > 1.0 && b > 1.0
        r = (pp_approx^2 - 3.0)/6.0
        s = 1.0/(2*aa - 1.0)
        t = 1.0/(2*bb - 1.0)
        h = 2.0/(s+t)
        w = pp_approx*sqrt(h+r)/h - (t-s)*(r + 5.0/6.0 - 2.0/(3.0*h))
        x = aa/ (aa+bb*exp(w^2))
    else
        r = 2.0*bb
        t = 1.0/(9.0*bb)
        t = r*(1.0-t+pp_approx*sqrt(t))^3
        if t <= 0.0
            x = -expm1((log((1.0-pp)*bb)+lb)/bb)
        else
            t = (4.0*aa+r-2.0)/t
            if t <= 1.0
                x = exp((log(pp*aa)+lb)/aa)
            else
                x = 1.0 - 2.0/(t+1.0)
            end
        end
    end

    #solve x using modified newton-raphson iteration

    r = 1.0 - aa
    t = 1.0 - bb
    pp_approx_prev = 0.0
    sq = 1.0
    prev = 1.0

    if x < 0.0001
        x = T(0.0001)
    end
    if x > .9999
        x = T(.9999)
    end

    iex = max(-5.0/aa^2 - 1.0/pp^0.2 - 13.0, -30.0)
    acu = 10.0^iex

    #iterate
    while true
        #pp_approx = beta_inc(aa,bb,x)[1]
        pp_approx = beta_inc_(aa,bb,x)
        xin = x
        pp_approx = (pp_approx-pp)*exp(lb+r*log(xin)+t*log1p(-xin))
        if pp_approx * pp_approx_prev <= 0.0
            prev = max(sq, fpu)
        end
        g = 1.0

        tx = x - g*pp_approx
        while true
            adj = g*pp_approx
            sq = adj^2
            tx = x - adj
            if (prev > sq && tx >= 0.0 && tx <= 1.0)
                break
            end
            g /= 3.0
        end

        #check if current estimate is acceptable

        if prev <= acu || pp_approx^2 <= acu
            x = tx
            return indx ? 1.0 - x : x
        end

        if tx == x
            return indx ? 1.0 - x : x
        end

        x = tx
        pp_approx_prev = pp_approx
    end
end



