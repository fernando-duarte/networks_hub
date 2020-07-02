using LinearAlgebra, Random, Distributions, NLsolve
using Test


function clearing_vec!(F, p, x)
    F .= p.-min.(p_bar,max.((1+g)*(A0'*p .+c .- x) .-g.*p_bar,0))
end

function shock()
    # draw shocks that produce the probability of default δ
    a = 1
    b = log.(δ)./(log.(1 .-w./c))
    dist = Beta.(a,[b...])
    draws = rand.(dist, 1)
    vcat(draws...).*c
    #rand(Float64,N) # normal
end

function loss(x)
     # total loss in actual network for shock x
    p = nlsolve((F, p)->clearing_vec!(F, p, x),[p_bar...], autodiff = :forward)
    @test converged(p)
    sum(x.+p_bar.-p.zero)
end

function loss_d(x)
    # total loss in disconnected network for shock x
    #b_d = b .+ max.(0, c.-b.-w ) # add fictitious outside liabilities
    #_d = c .+ max.(0, w.-c.+b ) # add fictitious outside assets
    sum(x+max.(0,x.-w))
end

## Glasserman and Young

y=0;
p_bar = (55.0+y, 55.0+y, 140.0, 55.0+y, 55.0+y); # total liabilities
c = (50.0, 50.0, 150.0, 50.0, 50.0); # outside assets
w = (5.0, 5.0, 10.0, 5.0, 5.0); # net worth
A0 = [0 y/p_bar[1] 0 0 0; 0 0 0 y/p_bar[1] 0; 10.0/p_bar[3] 10.0/p_bar[3] 0 10.0/p_bar[3] 10.0/p_bar[3]; 0 0 0 0 y/p_bar[1]; y/p_bar[1] 0 0 0 0]; # matrix of relative liabilities
g = 0.0; # bankruptcy costs
δ = (0.01, 0.1, 0.3, 0.2, 0.03); # probability of default
b = (55.0, 55.0, 100.0, 55.0, 55); # outside liabilities

a = w .+ p_bar; # total assets
d=  a .- c;# inside assets
f = p_bar .- b;# inside liabilities
β = (p_bar.-b)./p_bar # financial connectivity: proportion of liabilities inside the network
β⁺ = maximum(β)
N = length(c); # number of nodes

# example from Fernando's presentation and Glasserman-Young
x1 = [0.0, 0.0, 94.0, 0.0, 0.0]
@test loss(x1)==182.0
@test loss_d(x1)==178.0

ratio = loss(x1)/loss_d(x1)
bound = 1 + sum(δ.*c)/((1-β⁺)*sum(c))

@show ratio
@show bound

# draw shock from Beta distribution
x1 = shock()
ratio = loss(x1)/loss_d(x1)
bound = 1 + sum(δ.*c)/((1-β⁺)*sum(c));
@show ratio
@show bound

# example in https://github.com/siebenbrunner/NetworkValuation
c = (10.0,10.0,10.0); # net worth
matL = [0.0 20.0 10.0; 5.0 0.0 24.0; 10.0 0.0 0.0]
p_bar = [sum(matL,dims=2)...]
A0 = matL ./ repeat(p_bar, outer=[1,length(c)])
vecEquity = c .+ A0'*p_bar .- p_bar

clearing_vec!([0,0,0], p_bar, [0,0,0])
x0=[0.0,0.0,0.0]
p = nlsolve((F, p)->clearing_vec!(F, p, [x0...]),[p_bar...], autodiff = :forward)

true_p=[24.545454545454547,26.363636363636363,10.0]
@test norm(true_p-p.zero)<1e-10
