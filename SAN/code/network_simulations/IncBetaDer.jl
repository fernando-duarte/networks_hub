module IncBetaDer

# Goal: derivative of incomplete Beta function wrt the α and β parameters, and
# the same for the inverse of the incomplete Beta function. 
#
# I(x;α,β) with x ∈ [0, 1] gives the cdf of a Beta(α, β) distribution, whereas
# Iinv(x;α,β) with x ∈ [0, 1] gives the x-quantile for a Beta(α,β) distribution
# 
# I don't think there is an algorithm for the gradient of the inverse Beta
# function wrt α and β available, (note that rules for the gradients of inverse
# functions cannot be used for this). Of course, the numerical method for 
# the inverse Beta is a numerical inversion of the incomplete Beta function,
# which mainly calls `beta_inc`, so having an AD primitive for `beta_inc` is
# already nice.
#
# I will use notation from Boik and Robinson-Cox from here, i.e. α ≡ p, β ≡ q.

using SpecialFunctions
using ForwardDiff
import Base.Math: @horner
import SpecialFunctions: log1p
import ForwardDiff: Dual, value, partials

export beta_inc_grad, beta_inc_inv_, beta_inc_, discretize_beta

function discretize_beta(a, b, K)
    qstart = 1.0/(2.0 * K)
    qend = 1. - qstart
    xs = [beta_inc_inv_(a, b, x) for x in qstart:(1/K):qend]
    xs .* ((K*a/(a + b))/sum(xs))
end

include("beta_inc_grad.jl")
include("beta_inc_inv.jl")

end # module
