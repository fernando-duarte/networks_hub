module ExtendedPowerDist
# Implements the pdf, cdf, quantiles, and related functions, for the 
# Extended Power Distribution using MeasureTheory.jl
#= latex citation:
@misc{ogbonnaya2017extended,
      title={The extended power distribution: A new distribution on $(0, 1)$}, 
      author={Chibueze E. Ogbonnaya and Simon P. Preston and Andrew T. A. Wood},
      year={2017},
      eprint={1711.02774},
      archivePrefix={arXiv},
      primaryClass={stat.ME}
}
=#
export GenExtPow,ExtPow,LinFail,logdensity,cdf,logcdf,density,moment,quantile,mode

# Pkg.add(name="KeywordCalls", version="0.2.2") 
# Pkg.add(name="MeasureTheory", version="0.8.2") 
using MeasureTheory

# generalized extended power distribution with k parameters 
# α₀>0, α₁>=0, ..., αₖ>=0, and x in (0,1)
# case αᵢ==0 for i>=2 reduces to (non-generalized) extended power distribution ExtPow(α₀,α₁)
@parameterized GenExtPow(α) ≪ Lebesgue(𝕀)
function density(d::GenExtPow{(:α,)}, x)
    return (1/x)*sum((-1)^(h-1) * d.α[h] * h * log(x)^(h-1) for h=1:length(d.α))*cdf(GenExtPow(d.α),x) 
end

function logdensity(d::GenExtPow{(:α,)}, x)
    return -log(x)+log(sum((-1)^(h-1) * d.α[h] * h * log(x)^(h-1) for h=1:length(d.α))) + logcdf(GenExtPow(d.α),x) 
end

function cdf(d::GenExtPow{(:α,)}, x)
    # α₀==d.α[1], α₁==d.α[2], ..., d.αₖ==α[k+1]
    return exp(sum( (-1)^(h-1)*d.α[h]*log(x)^h for h=1:length(d.α) ))
end

function logcdf(d::GenExtPow{(:α,)}, x)
    # α₀==d.α[1], α₁==d.α[2], ..., d.αₖ==α[k+1]
    return sum((-1)^(h-1)*d.α[h]*log(x)^h for h=1:length(d.α))
end

# extended power distribution
# α₀>0, α₁>=0, and x in (0,1)
@parameterized ExtPow(α₀,α₁) ≪ Lebesgue(𝕀)
function density(d::ExtPow{(:α₀,:α₁)}, x)
    return ( d.α₀ - 2d.α₁*log(x) ) / x * exp( d.α₀ * log(x) - d.α₁ * log(x)^2)
end

function logdensity(d::ExtPow{(:α₀,:α₁)}, x)
    return  log( d.α₀ - 2d.α₁*log(x) )  + (d.α₀-1) * log(x) - d.α₁ * log(x)^2
end

function cdf(d::ExtPow{(:α₀,:α₁)}, x)
    return  exp( d.α₀ * log(x) - d.α₁ * log(x)^2)
end

function logcdf(d::ExtPow{(:α₀,:α₁)}, x)
    return  d.α₀ * log(x) - d.α₁ * log(x)^2
end

using SpecialFunctions # for erfc
function moment(d::ExtPow{(:α₀,:α₁)}, k)
    return  1 - (k/2)*sqrt(pi / d.α₁)*exp( (d.α₀ + k)^2 / (4d.α₁))*erfc( (d.α₀+k) / (2sqrt(d.α₁)) )
end

function quantile(d::ExtPow{(:α₀,:α₁)}, p)
    return   exp( (d.α₀ - sqrt(d.α₀^2 - 4d.α₁*log(p))) / (2d.α₁) )
end

function mode(d::ExtPow{(:α₀,:α₁)})
    return   ( (2d.α₀-1) - sqrt(1+ 8d.α₁) )  / (4d.α₁) 
end

# sampler 
# T ~ ExtPow(α₀,α₁) sample U~Uniform() and
# T = exp( (d.α₀ - sqrt(d.α₀^2 - 4d.α₁ log(U))) / (2d.α₁) ) if α₁ !=0
#     T = U^(1/d.α₀) otherwise 
# from cdf(ExtPow(α₀,α₁), T)==U, solve for T
        

# linear failure rate distribution
# α₀>0, α₁>=0, 0 < x < Inf
# If T~ExtPow(α₀,α₁) then -log(T) ~ LinFail(α₀,α₁)
@parameterized LinFail(α₀,α₁) ≪ Lebesgue(ℝ₊)
function density(d::LinFail{(:α₀,:α₁)}, x)
    return  (d.α₀ + 2d.α₁*x )*exp( -d.α₀*x - d.α₁*x^2)
end

end # module

#= tests
using ForwardDiff, Zygote, Quadrature, Cuba, Test, QuadGK

@testset "ExtPow" begin
    xgrid = 0.0001:0.01:(0.9999)
    for i=1:20
        α₀ = abs(randn())
        d1 = [density(ExtPow(α₀,0),x) for x in xgrid]
        d2 = [ MeasureTheory.density(Beta(α₀,1),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12

        d1 = [density(ExtPow(1,0),x) for x in xgrid]
        d2 = [MeasureTheory.density(Uniform(),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12

        d1 = [density(LinFail(α₀,0),x) for x in xgrid]
        d2 = [MeasureTheory.density(Exponential(α₀),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12

        α₁ = abs(randn())
        d1 = [density(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [density(GenExtPow([α₀,α₁]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        d1 = [density(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [density(GenExtPow([α₀,α₁,0]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        d1 = [density(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [density(GenExtPow([α₀,α₁,0,0,0]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        d1 = [cdf(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [cdf(GenExtPow([α₀,α₁]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        d1 = [cdf(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [cdf(GenExtPow([α₀,α₁,0]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        d1 = [cdf(ExtPow(α₀,α₁),x) for x in xgrid]
        d2 = [cdf(GenExtPow([α₀,α₁,0,0,0]),x) for x in xgrid]
        @test sum(abs2,d1-d2) ≈ 0 atol=1e-12
        
        sol = quadgk((x)->density(ExtPow(α₀,α₁),x),0,1; rtol=1e-18, atol=0, maxevals=10^7, order=7)
        @test sol[1] ≈ 1.0 atol=1e-8
        
        sol = quadgk((x)->density(LinFail(α₀,α₁),x),0,50; rtol=1e-18, atol=0, maxevals=10^7, order=7)
        @test sol[1] ≈ 1.0 atol=1e-8
       
        sol = quadgk((x)->density(GenExtPow([α₀,α₁]),x),0,1; rtol=1e-18, atol=0, maxevals=10^7, order=7)
        @test sol[1] ≈ 1.0 atol=1e-8
        
        sol = quadgk((x)->density(GenExtPow([α₀,α₁,0,0,0]),x),0,1; rtol=1e-18, atol=0, maxevals=10^7, order=7)
        @test sol[1] ≈ 1.0 atol=1e-8
        
        sol = quadgk((x)->density(GenExtPow([α₀,α₁,0.32,1.223,9.312]),x),0,1; rtol=1e-18, atol=0, maxevals=10^7, order=7)
        @test sol[1] ≈ 1.0 atol=1e-8
             
    end
end
=#

#= plots
using Plots: plot
xgrid = 0:0.01:1;
ygrid1 = [density(ExtPow(0.1,0.2),x) for x in xgrid];
ygrid2 = [density(ExtPow(1.1,0.2),x) for x in xgrid];
ygrid3 = [density(ExtPow(0.1,1.2),x) for x in xgrid];
ygrid4 = [density(ExtPow(1.1,1.2),x) for x in xgrid];
plot(xgrid, hcat(ygrid1,ygrid2,ygrid3,ygrid4))

xgrid = 0:0.01:1;
ygrid1 = [density(GenExtPow([0.1,0.2,0.1]),x) for x in xgrid];
ygrid2 = [density(GenExtPow([8,10.2,0.1]),x) for x in xgrid];
ygrid3 = [density(GenExtPow([0.1,1.2,5,20,20]),x) for x in xgrid];
ygrid4 = [density(GenExtPow([5.1,6.2,4.5,8.10]),x) for x in xgrid];
plot(xgrid, hcat(ygrid1,ygrid2,ygrid3,ygrid4))

xgrid = 0:0.05:10;
ygrid1 = [density(LinFail(0.1,0.2),x) for x in xgrid];
ygrid2 = [density(LinFail(1.1,0.2),x) for x in xgrid];
ygrid3 = [density(LinFail(0.1,1.2),x) for x in xgrid];
ygrid4 = [density(LinFail(1.1,1.2),x) for x in xgrid];
plot(xgrid, hcat(ygrid1,ygrid2,ygrid3,ygrid4))
=#

# compatible with symbolics
#= 
using Symbolics

@variables α,β,x
pdf_sym = density(ExtPow(α,β),x)
cdf_sym = cdf(ExtPow(α,β),x)
grad_pdf_sym = Symbolics.gradient(pdf_sym,vcat(α,β,x); simplify=true)
grad_cdf_sym = Symbolics.gradient(cdf_sym,vcat(α,β,x); simplify=true)
hess_pdf_sym = Symbolics.hessian(pdf_sym,vcat(α,β,x); simplify=true)
hess_cdf_sym = Symbolics.hessian(cdf_sym,vcat(α,β,x); simplify=true)
@test simplify(Symbolics.derivative(cdf_sym,x; simplify=true) - pdf_sym, expand=true)==0
=#

