import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()

using Convex, SCS, LinearAlgebra, Distributions, Test, Random, MathOptInterface
const MOI = MathOptInterface

# randomly generate indices for A where it is allowed to be non-zero
function randIndA(N::Integer)    
    lt = LowerTriangular(ones(N,N))-Matrix(I,N,N) 
    ind = findall( lt .> 0 )
    indᵀ = map(x->CartesianIndex(x[2],x[1]),CartesianIndices(lt)[ind])
    d = Distributions.Binomial(1,1/2)
    v = rand(d,length(ind))
    return [v[i]*LinearIndices(lt)[ind[i]] + (1-v[i])*LinearIndices(lt)[indᵀ[i]] for i=1: length(ind)]
end

function getA(N::Integer,nonzero_ind::Union{Nothing,Vector{Int}}=nothing)
    if isnothing(nonzero_ind)
        nonzero_ind = randIndA(N)
    end    
    zeros_ind = setdiff(1:N^2,nonzero_ind)
    A = Variable(N,N,Positive()) 
    problem = minimize(0.0)
    problem.constraints = [
                            A >= 0.0, A <= 1.0,
                            A[zeros_ind] == 0,
                            #b - p_bar + diagm(p_bar)*A*ones(N) == 0,
                            #c - a + transpose(A)*p_bar == 0
                           ] 
    Convex.solve!(problem,  SCS.Optimizer(verbose=false, eps=1e-12,linear_solver=SCS.DirectSolver), warmstart=false; silent_solver=true)
    if !(problem.status == MOI.INFEASIBLE || problem.status == MOI.ALMOST_INFEASIBLE)
        return evaluate(A)
    else
        return nothing
    end
end

function trainingSetA(N::Integer; sample_size::Integer = 20, max_iter::Integer = 2000)
    q = 0 # number of matrices found after max_iter iterations
    sampleA = Array{Array{Float64,2},1}(undef,sample_size) # initialize
    iter = 1 
    while iter<=max_iter && q<sample_size
        newA = getA(N)
        if !isnothing(newA)
            sampleA[q+1] = newA
            q+=1
        end
        iter+=1
    end
    return sampleA[1:q]
end

# test on network from Glasserman-Young paper
b = [55.0,  55.0,  100.0,  55.0,  55.0] # outside liabilities
c = [50.0, 50.0, 150.0, 50.0, 50.0] # outside assets
w = [5.0, 5.0, 10.0, 5.0, 5.0]; # net worth
p_bar = [55.0, 55.0, 140.0, 55.0, 55.0]; # total liabilities
a = w .+ p_bar; # total assets
N=length(b)

trueA = [0 0 0 0 0; 0 0 0 0 0; 10.0/p_bar[3] 10.0/p_bar[3] 0 10.0/p_bar[3] 10.0/p_bar[3]; 0 0 0 0 0; 0 0 0 0 0]
nonzero_ind = LinearIndices(trueA)[findall(trueA.>0)]
@test all(trueA .>= 0)
@test all(trueA .<= 1)
@test all([trueA[i,j]*trueA[j,i] for i=1:N for j=1:i] .== 0)
# @test all(b - p_bar + diagm(p_bar)*trueA*ones(N).== 0)
# @test all(c - a + transpose(trueA)*p_bar.== 0)
# solA = getA(5,nonzero_ind)
# @test solA ≈ trueA

# generate a number Q of matrices A that are consistent with b, c, w, p_bar, a
max_iter = 10000
Q = 20 # number of matrices we want to generate
A_in = trainingSetA(N,sample_size=Q,max_iter=max_iter)

# check generated matrices are consistent with b, c, w, p_bar, a
q=length(A_in)
if q>0
    ε = 1e-8
    @testset "checkA" begin
        for r=1:q
            all(A_in[r] .>= 0-ε)
            all(A_in[r] .<= 1+ε)
            [ (@test A_in[r][i,j]*A_in[r][j,i] ≈ 0 atol=ε) for i=1:N for j=1:i ]
            #[ (@test (b - p_bar + diagm(p_bar)*A_in[r]*ones(N))[i] ≈ 0 atol=ε) for i=1:N ]
            #[ (@test (c - a + transpose(A_in[r])*p_bar)[i] ≈ 0 atol=ε) for i=1:N ]
        end
    end
end

# test on empirical network

# N = 5
# b= [0.8774984375, 0.9302580625, 0.8084824375, 0.7275786875, 0.721182]
# c=[1.451636125, 1.38565475, 1.516806125, 1.158976375, 0.721182]
# w=[0.168059, 0.144022, 0.17871046875, 0.102316, 0.04922]
# p_bar= [2.006993, 1.794448, 1.6433575, 1.207323, 1.032553]
# a=[2.175052, 1.93847, 1.82206796875, 1.309639, 1.081773]
# N=length(b)

# N=20
b= [0.8774984375, 0.9302580625, 0.8084824375, 0.7275786875, 0.7211820000000001, 0.6814913333333333, 0.6582686666666667, 0.5838533333333333, 0.232036765625, 0.306932, 0.20765666666666668, 0.177556671875, 0.1259010859375, 0.1000723515625, 0.110032796875, 0.0901028203125, 0.115514, 0.0946800625, 0.10999406266666667, 0.0826875234375]
c=[1.451636125, 1.38565475, 1.516806125, 1.158976375, 0.7211820000000001, 0.6814913333333333, 0.6582686666666667, 0.5838533333333333, 0.37314634375, 0.306932, 0.20765666666666668, 0.247988171875, 0.249620546875, 0.128101765625, 0.171698203125, 0.0915575625, 0.115514, 0.14428496875, 0.10999406266666667, 0.143162515625]
w=[0.168059, 0.144022, 0.17871046875, 0.102316, 0.04922, 0.082995, 0.035765, 0.0406288125, 0.02398659375, 0.018699, 0.012557, 0.027647849609375, 0.027033, 0.028054, 0.02248987109375, 0.012782736328125, 0.009780671875, 0.02661243359375, 0.005355671875, 0.0160814091796875]
p_bar=  [2.006993, 1.794448, 1.6433575, 1.207323, 1.032553, 0.939242, 0.951638, 0.835151188, 0.47769165625, 0.441699, 0.298928, 0.26344503125, 0.239999, 0.209598, 0.16664809375, 0.16384959375, 0.16349032800000002, 0.139301015625, 0.159635422, 0.135933609375]
a=[2.175052, 1.9384700000000001, 1.82206796875, 1.309639, 1.081773, 1.022237, 0.987403, 0.8757800005, 0.50167825, 0.46039800000000003, 0.311485, 0.291092880859375, 0.267032, 0.237652, 0.18913796484375, 0.17663233007812498, 0.17327099987500003, 0.16591344921875, 0.164991093875, 0.1520150185546875]
N=length(b)

max_iter = 10000
Q = 20 # number of matrices we want to generate
A_in = trainingSetA(N,sample_size=Q,max_iter=max_iter)

# check generated matrices are consistent with b, c, w, p_bar, a
q=length(A_in)
if q>0
    ε = 1e-8
    @testset "checkA" begin
        for r=1:q
            all(A_in[r] .>= 0-ε)
            all(A_in[r] .<= 1+ε)
            [ (@test A_in[r][i,j]*A_in[r][j,i] ≈ 0 atol=ε) for i=1:N for j=1:i ]
            #[ (@test (b - p_bar + diagm(p_bar)*A_in[r]*ones(N))[i] ≈ 0 atol=ε) for i=1:N ]
            #[ (@test (c - a + transpose(A_in[r])*p_bar)[i] ≈ 0 atol=ε) for i=1:N ]
        end
    end
end