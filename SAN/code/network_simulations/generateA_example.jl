import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()

using Convex, SCS, LinearAlgebra, Distributions, Test, Random, MathOptInterface, Parameters
const MOI = MathOptInterface
include("generateA.jl")
using .generateA

## create a sample of Q matrices A that satisfy A[i,j]*A[j,i]==0 and 0 .<= A .<= 1
max_iter = 10000
Q = 20 
A_in_net5 = trainingSetA(5,sample_size=Q,max_iter=max_iter,con="netting") # N=5 nodes, can use for training NN
A_in_net20 = trainingSetA(20,sample_size=Q,max_iter=max_iter,con="netting") # N=20 nodes, can use for training NN

# can also use "perturbed" matrices
δ = 1e-4
R = 10 # number of perturbed matrices for each matrix in A_in
Aδ_in_net5 = Matrix{Matrix{Float64}}(undef,size(A_in_net5,1),R)
for q=1:Q
    for r=1:R
        Aδ_in_net5[q,r] = min.(1, max.(0, A_in_net5[q] + rand(5,5)*δ))      
    end
end
Aδ_in_net5 = reshape(Aδ_in_net5,Q*R) # can be used for training NN

#=
#####################################################################
Testing and experiments
#####################################################################
=#

# test that generated matrices satisfy A[i,j]*A[j,i]==0 and 0 .<= A .<= 1
if length(A_in_net5)>0
    ε = 1e-12
    A_in = A_in_net5
    N = 5
    q=length(A_in)
    @testset "check_net5" begin
        for r=1:q
            all(A_in[r] .>= 0-ε)
            all(A_in[r] .<= 1+ε)
            [ (@test A_in[r][i,j]*A_in[r][j,i] ≈ 0 atol=ε) for i=1:N for j=1:i ]
            #[ (@test (b - p_bar + diagm(p_bar)*A_in[r]*ones(N))[i] ≈ 0 atol=ε) for i=1:N ]
            #[ (@test (c - a + transpose(A_in[r])*p_bar)[i] ≈ 0 atol=ε) for i=1:N ]
        end
    end
end

if length(A_in_net20)>0
    ε = 1e-12
    A_in = A_in_net20
    N = 20
    q=length(A_in)
    @testset "check_net20" begin
        for r=1:q
            all(A_in[r] .>= 0-ε)
            all(A_in[r] .<= 1+ε)
            [ (@test A_in[r][i,j]*A_in[r][j,i] ≈ 0 atol=ε) for i=1:N for j=1:i ]
            #[ (@test (b - p_bar + diagm(p_bar)*A_in[r]*ones(N))[i] ≈ 0 atol=ε) for i=1:N ]
            #[ (@test (c - a + transpose(A_in[r])*p_bar)[i] ≈ 0 atol=ε) for i=1:N ]
        end
    end
end

## try to find A that satisfy all constraints
# create net with N=5 from the data
@with_kw mutable struct Net5
    b= [0.8774984375, 0.9302580625, 0.8084824375, 0.7275786875, 0.721182]
    c= [1.451636125, 1.38565475, 1.516806125, 1.158976375, 0.721182]
    w= [0.168059, 0.144022, 0.17871046875, 0.102316, 0.04922]
    p_bar= [2.006993, 1.794448, 1.6433575, 1.207323, 1.032553]
    a= [2.175052, 1.93847, 1.82206796875, 1.309639, 1.081773]
    N= length(b)
end
net5 = Net5()

# create A with all constraints
max_iter = 10000
Q = 20 
A_in_net5 = trainingSetA(net5,sample_size=Q,max_iter=max_iter)
isnothing(A_in_net5) # didn't find any A

# test all possible zero/nonzero patterns for A
aI = allIndA(net5.N)
allA = Vector{Union{Matrix{Float64}, Nothing}}(undef, length(aI[1]))
allA = getA(net5,nonzero_ind=aI[1])
allA2 = filter(x->!isnothing(x),allA)
isempty(allA2) # there is no feasible A

## test on network from Glasserman-Young paper, where we know we can find A that satisfies all constraints
# create net
@with_kw mutable struct NetGY
    b = [55.0,  55.0,  100.0,  55.0,  55.0] # outside liabilities
    c = [50.0, 50.0, 150.0, 50.0, 50.0] # outside assets
    w = [5.0, 5.0, 10.0, 5.0, 5.0]; # net worth
    p_bar = [55.0, 55.0, 140.0, 55.0, 55.0]; # total liabilities
    a = w .+ p_bar; # total assets
    d = a .- c; # inside assets
    f = p_bar .- b; # inside liabilities
    N::Integer = length(b); # number of nodes
end
net = NetGY()

# generate a number Q of matrices A that are consistent with net
max_iter = 10000
Q = 20 # number of matrices we want to generate
A_in = trainingSetA(net,sample_size=Q,max_iter=max_iter) 
isnothing(A_in_net5) # found some A's
@show round.(A_in[1],digits=6) # display first matrix A found

# check generated matrices are consistent with net
q=length(A_in)
if q>0
    ε = 1e-8
    @testset "checkA" begin
        @unpack b,c,p_bar,a,N = net 
        for r=1:q
            all(A_in[r] .>= 0-ε)
            all(A_in[r] .<= 1+ε)
            [ (@test A_in[r][i,j]*A_in[r][j,i] ≈ 0 atol=ε) for i=1:N for j=1:i ]
            [ (@test (b - p_bar + diagm(p_bar)*A_in[r]*ones(N))[i] ≈ 0 atol=ε) for i=1:N ]
            [ (@test (c - a + transpose(A_in[r])*p_bar)[i] ≈ 0 atol=ε) for i=1:N ]
        end
    end
end

## let's compare with the A matrix from Glasserman and Young
gyA = [0 0 0 0 0; 0 0 0 0 0; 10.0/net.p_bar[3] 10.0/net.p_bar[3] 0 10.0/net.p_bar[3] 10.0/net.p_bar[3]; 0 0 0 0 0; 0 0 0 0 0] # the A matrix from GY
nonzero_indGY = LinearIndices(gyA)[findall(gyA.>0)] # indices of non-zero elements

# check gyA satisfies constraints
@testset "gyA" begin
    @unpack b,c,p_bar,a,N = net 
    @test all(gyA .>= 0)
    @test all(gyA .<= 1)
    @test all([gyA[i,j]*gyA[j,i] for i=1:N for j=1:i] .== 0)
    @test all(b - p_bar + diagm(p_bar)*gyA*ones(N).== 0)
    @test all(c - a + transpose(gyA)*p_bar.== 0)
end

# use `getA` to create a Matrix with the same pattern of non-zero elements
solA = getA(net,nonzero_ind=nonzero_indGY)
@test solA ≈ gyA

# test all possible zero/nonzero patterns for A
aI = allIndA(net.N)
allA = Vector{Union{Matrix{Float64}, Nothing}}(undef, length(aI[1]))
allA = getA(net,nonzero_ind=aI[1])
allA2 = filter(x->!isnothing(x),allA)
all([allA2[i] ≈ gyA for i=1:length(allA2)]) # the only feasible matrix is the gyA matrix
