include("NU2_preamble.jl")

f = eval(NonLinProbPrecompile.f_noeval_good[1])
f_num = eval(NonLinProbPrecompileNum.f_noeval_num[1])
f_obj = eval(NonLinProbPrecompileObj.f_noeval_obj)
@variables z[1:M]
p=[]
obj_expr = f_obj(z,p)
f_expr = f_num(z,p)
Q = length(f_expr) # number of constraints
eq = 0 .~ f_expr
ns = NonlinearSystem(eq,z,p)
os = OptimizationSystem(obj_expr,z,[])

obj_num = generate_function(os)
grad_num = generate_gradient(os)
sys_num = generate_function(ns) 
jac_num = generate_jacobian(ns, sparse=true)
jac_sym = generate_jacobian(ns,z,p)
hess_obj_num = generate_hessian(os, sparse=true)
hess_obj_sym = calculate_hessian(os)
hess_sym(i) = Symbolics.sparsehessian(f_expr[i],Symbolics.scalarize(z))

on = @eval eval(obj_num)
gn = @eval eval(grad_num[1])
gn_iip = @eval eval(grad_num[2])
sn = @eval eval(sys_num[1])
sn_iip = @eval eval(sys_num[2])
jn = @eval eval(jac_num[1])
js_iip = @eval eval(jac_sym[2])
p0=Float64[];
@show on(z0,p0)
@show gn(z0,p0)
@show sn(z0,p0)
@show jn(z0,p0)

gn_z0 = zeros(length(z0)) 
sn_z0 = zeros(length(eq)) 
jn_z0 = spzeros(length(eq),length(z0)) 
@show gn_iip(gn_z0,z0,p0)
@show sn_iip(sn_z0,z0,p0)
@show js_iip(jn_z0,z0,p0)

@show gn_z0
@show sn_z0
@show jn_z0

@show jac_sp = ModelingToolkit.jacobian_sparsity(ns);
@show jac_sp.colptr;
@show jac_sp.rowval;

jac_nzval = calculate_jacobian(ns)[jac_sp]
jac_nzval_num = build_function(jac_nzval,z) 
jac_nzval_oop = @eval eval(jac_nzval_num[1])
jac_nzval_iip = @eval eval(jac_nzval_num[2])
jac_nzval0 = spzeros(nnz(jac_sp)) 
jac_nzval_oop(z0)
jac_nzval_iip(jac_nzval0,z0)
@show jac_nzval0

rows_jac_sp = jac_sp.rowval
cols_jac_sp = vcat([fill(j,length(nzrange(jac_sp,j))) for j=1:M]...)
@test sparse(rows_jac_sp,cols_jac_sp,true,Q,M)==jac_sp

# sparsity pattern of the Lagrangian's hessian
for i=1:Q
   @eval @parameters ($(Symbol("lambda$i"))) 
end
λvec = [ @eval ($(Symbol("lambda$i"))) for i=1:Q]
@parameters σ

lag =  σ * obj_expr + dot(f_expr,λvec) 
temp_sp = OptimizationSystem(lag,z,vcat(λvec...,σ))
hess_sp = ModelingToolkit.hessian_sparsity(temp_sp)
rows_hess_sp_full = hess_sp.rowval
cols_hess_sp_full = vcat([fill(j,length(nzrange(hess_sp,j))) for j=1:M]...)
@test sparse(rows_hess_sp_full,cols_hess_sp_full,true,M,M)==hess_sp

hess_sp_lt = sparse(LowerTriangular(Array(hess_sp)))
rows_hess_sp = hess_sp_lt.rowval
cols_hess_sp = vcat([fill(j,length(nzrange(hess_sp_lt,j))) for j=1:M]...)

@parameters λ[1:Q] σ
λcat = vcat(λ...)
hess_λ = σ * hess_obj_sym + sum(λ[i]*hess_sym(i) for i=1:Q) 
hess_num = build_function(hess_λ,z,λ,σ) # function (z, λ, σ) or function (var,z, λ, σ)
hess_oop = @eval eval(hess_num[1])
hess_iip = @eval eval(hess_num[2])
hess_z0 = spzeros(length(z0),length(z0)) 
λ0 = ones(Q)
σ0 = 1.0
hess_oop(z0,λ0,σ0)
hess_iip(hess_z0,z0,λ0,σ0)

hess_σ = σ * hess_obj_sym
hess_σ_num = build_function(hess_σ,z,σ)
hess_σ_oop = @eval eval(hess_σ_num[1])
hess_σ_iip = @eval eval(hess_σ_num[2])
hess_σ_z0 = spzeros(length(z0),length(z0)) 
σ0 = 1.0
hess_σ_oop(z0,σ0)
hess_σ_iip(hess_σ_z0,z0,σ0)
@show hess_σ_z0

hess_σ_nzval = hess_σ[hess_sp_lt]
hess_σ_nzval_num = build_function(hess_σ_nzval,z,σ) 
hess_σ_nzval_oop = @eval eval(hess_σ_nzval_num[1])
hess_σ_nzval_iip = @eval eval(hess_σ_nzval_num[2])
hess_σ_nzval0 = spzeros(nnz(hess_sp_lt)) 
σ0 = 1.0
hess_σ_nzval_oop(z0,σ0)
hess_σ_nzval_iip(hess_σ_nzval0,z0,σ0)
@show hess_σ_nzval0

hess_λ_nzval = hess_λ[hess_sp_lt]
hess_λ_nzval_num = build_function(hess_λ_nzval,z,λ,σ) # function (z, λ, σ) or function (var,z, λ, σ)
hess_λ_nzval_oop = @eval eval(hess_λ_nzval_num[1])
hess_λ_nzval_iip = @eval eval(hess_λ_nzval_num[2])
hess_λ_nzval0 = spzeros(nnz(hess_sp_lt)) 
λ0 = ones(Q)
σ0 = 1.0
hess_λ_nzval_oop(z0,λ0,σ0)
hess_λ_nzval_iip(hess_λ_nzval0,z0,λ0,σ0)
@show hess_σ_nzval0

@parameters v[1:M]
v_cat = vcat(v...)
Hv_sym = hess_λ*v_cat
Hv_num = build_function(Hv_sym,z,λ,v,σ)
Hv_oop = @eval eval(Hv_num[1])
Hv_iip = @eval eval(Hv_num[2])
Hv0 = spzeros(M) 
v0 = ones(M)
Hv_oop(z0,λ0,v0,σ0)
Hv_iip(Hv0,z0,λ0,v0,σ0)
@show Hv0
@test Array(Hv0)==Hv_oop(z0,λ0,v0,σ0)

############################
include("Network_NLPModel.jl")
using  NLPModelsIpopt

nlp = NetNLP(z0,low,upp,Q,jac_sp,hess_sp,N; minimize = true,name = "Network Optimization",p=[])
stats = ipopt(nlp,output_file="ipopt_out.txt")

                           
################ tests ###############
using ForwardDiff, Zygote, FiniteDiff
c_vals = zeros(Q)
hess_vals = zeros(nnz(hess_sp_lt))
jac_vals = zeros(nnz(jac_sp));jac_rows=zeros(Int,nnz(jac_sp));jac_cols=similar(jac_rows)
hess_vals = zeros(nnz(hess_sp_lt));hess_rows=zeros(Int,nnz(hess_sp_lt));hess_cols=similar(hess_rows)

function __c(z,p,x)
    vcat(
        c_lin(z,p),
        c_quad(z,p),
        c_p(z,p,x)
        )
end

p0 = p_dense
@test NetDefs.obj(z0,p0)==NLPModels.obj(nlp,z0)
@test vcat(NetDefs.c_lin(z0,p0), NetDefs.c_quad(z0,p0), NetDefs.c_p(z0,p0,x0)) == __c(z0,p0,x0)
@test sum(abs2,__c(z0,p0,x0)-NLPModels.cons!(nlp,z0,c_vals)) ≈ 0 atol=1e-30
NLPModels.cons!(nlp,z0,c_vals);
@test sum(abs2,__c(z0,p0,x0)-c_vals) ≈ 0 atol=1e-30

rr, cc = NLPModels.jac_structure!(nlp,jac_rows,jac_cols)
vv = NLPModels.jac_coord!(nlp,z0, jac_vals)
jac1 = Array(sparse(rr,cc,vv))
jac2 = Array(sparse(jac_rows,jac_cols,jac_vals))
jac3 = FiniteDiff.finite_difference_jacobian(z->__c(z,p0,x0),z0)
@test jac1 ≈ jac2 ≈ jac3 


rr_h, cc_h = NLPModels.hess_structure!(nlp,hess_rows,hess_cols)

hess_σ_iip(vals,z,obj_weight)
vv_h = NLPModels.hess_coord!(nlp,z0,hess_vals;obj_weight=1)
hess1 = Array(sparse(rr_h,cc_h,vv_h))
hess2 = Array(sparse(hess_rows,hess_cols,hess_vals))
hess3 = FiniteDiff.finite_difference_hessian(z->NLPModels.obj(nlp,z),z0)
hess4 = FiniteDiff.finite_difference_hessian(z->NetDefs.obj(z,p0),z0)
@test hess1 ≈ hess2 ≈ hess3 ≈ hess4



pnew = varmap_to_vars([β=>3.0, c=>10.0, γ=>2.0],parameters(sys))






fun(z) =  0.0
function fun_grad!(g, x)
    g .= zeros(length(x))
end
function fun_hess!(h, x)
    h .= zeros(length(x),length(x))
end
df = TwiceDifferentiable(fun, fun_grad!, fun_hess!, z0)

con_c!(c, x) = sn_iip(c,x,[])
function con_jacobian!(J, x)
    js_iip(J, x,[])
    return J
end


function con_h!(h, x, λ)
   sn_iip(h, x, λ)
   return h
end

z00= low + (upp-low)/2
df = TwiceDifferentiable(fun, fun_grad!, fun_hess!, z00)
lc = spzeros(Q)
uc = spzeros(Q)
dfc = TwiceDifferentiableConstraints(con_c!, con_jacobian!, con_h!,low,upp, lc, uc)

res = optimize(df, dfc, z00, IPNewton(show_linesearch = true), 
    Optim.Options(
                    g_tol = 1e-12,
                    x_tol = 1e-12,
                    f_tol = 1e-12,
                    allow_f_increases = true,
                    successive_f_tol = 5,
                    show_trace = true,
                    iterations=100
                  )
                )

zsol = res.minimizer;
c_quad(zsol,p_dense)
c_lin(zsol,p_dense)
c_p(zsol,p_dense,x0)


z00=low + (upp-low)/2
sol_NM = optimize(fun, z0, NelderMead())
sol_LBFGS = optimize(fun, z0, LBFGS())
sol_LBFGSad =optimize(fun, z0, LBFGS(); autodiff = :forward)
sol_SAMIN =optimize(fun, low,upp,z00, SAMIN())


df = TwiceDifferentiable(fun, fun_grad!, fun_hess!, z00)
dfc = TwiceDifferentiableConstraints(low, upp)

res = optimize(df, dfc, z00, IPNewton(), 
    Optim.Options(
                    g_tol = 1e-12
                  )
                )


#######################NLPModelsIpopt.jl###################################

using ADNLPModels, NLPModels, NLPModelsIpopt


nlp = ADNLPModel(fun,z00,low,upp,con_c!, lc, uc)
fx = obj(nlp, nlp.meta.x0)
gx = grad(nlp, nlp.meta.x0)
Hx = hess(nlp, nlp.meta.x0)

##########################################################


include("NU2_preamble.jl")

f = eval(NonLinProbPrecompile.f_noeval_good[1])
f_num = eval(NonLinProbPrecompileNum.f_noeval_num[1])
@variables z[1:M]
# @parameters p[1:P]
# f_expr = f(z,p)
# eq = 0 .~ f(z,p)
p=[]
f_expr = f_num(z,p)
eq = 0 .~ f_expr
#p0 = p_dense
# ns = NonlinearSystem(eq,z,p; name=:random_name, defaults=merge(Dict(vcat(z) .=> z0),Dict(vcat(p) .=> p0)) ) 
ns = NonlinearSystem(eq,z,p)
# prob_ns = NonlinearProblem(ns,z0,p0; check_length=false,jac = true, sparse=true,checkbounds = false, linenumbers = false) 
#prob_ns = NonlinearProblem(ns,z0,p0; check_length=false)
# NonlinearFunction(ns,z0,p0)
@show states(ns);
@show parameters(ns);


sys_num = generate_function(ns) 
sys_sym = generate_function(ns,z,p) 
jac_num = generate_jacobian(ns, sparse=true, skipzeros=false)
jac_sym = generate_jacobian(ns,z,p) 
jac_sym_sp = Symbolics.sparsejacobian(f_expr,states(ns))
hess_sym = ModelingToolkit.hessian(f_expr[1],states(ns); simplify=true) # hessian of 1st equation in f_expr
hess_num = build_function(hess_sym, states(ns); expression = false, target = Symbolics.JuliaTarget())

Q = 2N+1+N*mini_batch #number of constraints
@variables λ[1:Q]
[λ[i]*ModelingToolkit.hessian(f_expr[i],states(ns); simplify=true) for i=1:Q]

hv_sym = 

sn = @eval eval(sys_num[1])
sn_iip = @eval eval(sys_num[2])
ss = @eval eval(sys_sym[1])

jn = @eval eval(jac_num[1])
jn_iip = @eval eval(jac_num[2])
js = @eval eval(jac_sym[1])
js_iip = @eval eval(jac_sym[2])
js_sp = jac_sym_sp

hn = @eval eval(hess_num[1])
hn_iip = @eval eval(hess_num[2])
hs = @eval eval(hess_sym)


p0=Float64[];
@show sn(z0,p0)
@show ss(z,p)
@show jn(z0,p0)
@show js(z,p)
@show js_sp
@show hn(z0,p0)
@show hs

sn_z0 = zeros(length(eq)) 
jn_z0 = spzeros(length(eq),length(z0)) 
hn_z0 = spzeros(length(z0),length(z0)) 
@show sn_iip(sn_z0,z0,p0)
@show jn_iip(jn_z0,z0,p0)
@show hn_iip(hn_z0,z0,p0)

@show sn_z0
@show jn_z0
@show hn_z0


fun(z) =  0.0
function fun_grad!(g, x)
    g .= zeros(length(x))
end
function fun_hess!(h, x)
    h .= zeros(length(x),length(x))
end
df = TwiceDifferentiable(fun, fun_grad!, fun_hess!, z0)

con_c!(c, x) = sn_iip(c,x,[])
function con_jacobian!(J, x)
    jn_iip(J, x,[])
    return J
end
function con_h!(h, x, λ)
    h[1,1] += λ[1]*2
    h[2,2] += λ[1]*2
end;



f_optim(z) = sn(z,p0)
optimize(f_optim, z0, NelderMead())
optimize(f_optim, z0, LBFGS())
optimize(f_optim, z0, LBFGS(); autodiff = :forward)
g!(G, x) = 


inner_optimizer = GradientDescent()
results = optimize((x)->sn(x,p0), low, upp, (upp-low)/2+low, Fminbox(inner_optimizer), autodiff = :forward)

nlsolve((F,x)->sn_iip(F,x,p0), z0, autodiff = :forward)
nlsolve((F,x)->sn_iip(F,x,p0), (F,x)->jn_iip(F,x,p0), z0)


######################################
c_lin(z0,p0)
c_quad(z0,p0)
c_chance(z0,p0)
_c(z0,p0)

c_p(z0,p0,x0)
_c(z0,p0,x0)

contraction(z0,p0,x0)
p1=p_contraction(z0, p0, x0, 10)
p2=p_nlsolve(z0,p0,x0)
p1 ≈ p2


obj(z0,p0)
obj_contraction(z0,p0,x0)
obj_nlsolve(z0,p0,x0)
α,β= selα*z0,selβ*z0
betaprod_pdf(x0,α,β)
betaprod_density(x0,α,β)
dx0=zeros(1,mini_batch)
int_nlsolve!(dx0,x0,vcat(z0,p0))
int_cont!(dx0,x0,vcat(z0,p0))
int_nlsolve(x0,vcat(z0,p0))
int_cont(x0,vcat(z0,p0))
obj_nlsolve_int(z0,p0)
obj_cont_int(z0,p0)
obj_nlsolve_intH(z0,p0)
obj_cont_intH(z0,p0)

obj_full(obj(z0,p0),p0)
obj_full(obj_contraction(z0,p0,x0),p0)
obj_full(obj_nlsolve(z0,p0,x0),p0)

# integrate by drawing x
dist_uni = Distributions.Product(Distributions.Uniform.(zeros(size(α)),ones(size(β))))
dist_beta = Distributions.Product(Distributions.Beta.(α,β))
num_draw = 100000
num_draw = (num_draw ÷ mini_batch)*mini_batch # a multiple of mini_batch

uni_draw = rand(dist_uni,num_draw)
beta_draw = rand(dist_beta,num_draw)
@show isapprox(mean(uni_draw,dims=2),fill(0.5,N), atol=1e-3);
@show isapprox(mean(beta_draw,dims=2),α./(α + β), atol=1e-3);

uni_block = BlockArray(uni_draw,[N],fill(mini_batch,num_draw ÷ mini_batch))
beta_block = BlockArray(beta_draw,[N],fill(mini_batch,num_draw ÷ mini_batch))

# using contraction mapping
int_cont0(x) = int_cont!(dx0,x,vcat(z0,p0))
obj_cont0(x) = obj_contraction(z0,p0,x)
approx_uni  = mean([mean(int_cont0(view(uni_block, Block(i)))) for i=1:blocklength(uni_block)])
approx_beta = mean([mean(obj_cont0(view(beta_block, Block(i)))) for i=1:blocklength(beta_block)])
approx_int  = obj_cont_int(z0,p0)
@show approx_uni, approx_beta, approx_int ;

# using nlsolve
int_nlsolve0(x) = int_nlsolve!(dx0,x,vcat(z0,p0))
obj_nlsolve0(x) = obj_nlsolve(z0,p0,x)
approx_uni = mean([mean(int_nlsolve0(view(uni_block, Block(i)))) for i=1:blocklength(uni_block)])
approx_beta =mean([mean(obj_nlsolve0(view(beta_block, Block(i)))) for i=1:blocklength(beta_block)])
approx_int = obj_nlsolve_int(z0,p0)
@show approx_uni, approx_beta, approx_int ;

# time cuba vs HCubature
using Test, BenchmarkTools
@btime obj_cont_int($z0,$p0)
@btime obj_cont_intH($z0,$p0)

# gradients and hessians of integrals
using Zygote, ForwardDiff, FiniteDiff 
testf(p) = obj_cont_int(z0,p)
# gradients
Zygote.gradient(testf,p0)
ForwardDiff.gradient(testf,p0)
FiniteDiff.finite_difference_gradient(testf,p0)

# hessians
@test_broken Zygote.hessian(testf,p0)
@test_broken ForwardDiff.hessian(testf,p0)
@test_broken Zygote.jacobian(p->Zygote.gradient(testf,p),p0)
@test_broken ForwardDiff.jacobian(x->ForwardDiff.gradient(testf,x),p0)
@test_broken Zygote.jacobian(p->ForwardDiff.gradient(testf,p),p0)
@test_broken ForwardDiff.jacobian(x->Zygote.gradient(testf,x),p0)
@test_broken Zygote.jacobian(x->FiniteDiff.finite_difference_gradient(testf,x),p0)

# these do work but can be slow 
@test_skip FiniteDiff.finite_difference_hessian(testf,p0)   
@test_skip ForwardDiff.jacobian(x->FiniteDiff.finite_difference_gradient(testf,x),p0)
@test_skip FiniteDiff.finite_difference_jacobian(x->ForwardDiff.gradient(testf,x),p0)
@test_skip FiniteDiff.finite_difference_jacobian(x->Zygote.gradient(testf,x)[1],p0)

## modling toolkit
using SparsityDetection, SparseArrays
input = rand(10)
output = similar(input)
sparsity_pattern = jacobian_sparsity(f,output,input)
jac = Float64.(sparse(sparsity_pattern))
using SparseDiffTools
colors = matrix_colors(jac)
using FiniteDiff
FiniteDiff.finite_difference_jacobian!(jac, f, rand(30), colorvec=colors)
@show fcalls # 5
forwarddiff_color_jacobian!(jac, f, x, colorvec = colors)



## Modeling Toolkit      #https://github.com/JuliaStats/StatsFuns.jl/tree/an/nopdf

module NonLinProbPrecompile
    using ModelingToolkit, LinearAlgebra

    function system(; kwargs...)
        # Define some variables
        M = 50;P=3031;x0=ones(M)/10;
    
        @variables z[1:M]
        @parameters p[1:P]
        zcat = vcat(z...);
        pcat = vcat(p...);
        # Define a system of nonlinear equations
        ceq = vcat(
            0 .~ c_lin(zcat,pcat),
            0  ~ c_quad(zcat,pcat) ,
            0 .~ c_p(zcat,pcat,x0)
          )
        ns = NonlinearSystem(ceq,z,p)
        return generate_function(ns,z,p)
    end
    # Setting eval_expression=false and eval_module=[this module] will ensure
    # the RGFs are put into our own cache, initialised below.
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
    const f_noeval_good = system(; eval_expression=false, eval_module=@__MODULE__)
end

f = eval(Main.NonLinProbPrecompile.f_noeval_good[1])
















function c_quad(z::Symbolics.Arr{Num},p::Symbolics.Arr{Num})
    H  = reshape(p[6N+2+2N*M:6N+1+2N*M+M^2],M,M)
    return dot(z,Symbolics.scalarize(H*z)) #transpose(z)*H*z
end

@variables z[1:M]
@parameters p[1:P]
zcat = vcat(z...);
pcat = vcat(p...);

ceq = vcat(
        0 .~ c_lin(zcat,pcat),
        0  ~ c_quad(zcat,pcat) ,
        0 .~ c_p(zcat,pcat,x0)
      )

loss = Symbolics.scalarize(Symbolics.scalarize(sum(selp*zcat)) + Symbolics.scalarize(sum(Symbolics.scalarize(selc*zcat).*x0)) )
sys = OptimizationSystem(loss,z,p; equality_constraints=ceq)
u0 = zcat .=> z0
p0 = pcat .=> p_dense
prob = OptimizationProblem(sys,u0,p0,grad=true,hess=true)  #
sol = solve(prob,IPNewton())

#@test ModelingToolkit.varmap_to_vars(u0,states(ns); defaults=ModelingToolkit.defaults(ns))==z0

ns = NonlinearSystem(ceq,z,p)
prob_ns = NonlinearProblem(ns,z0,p_dense; check_length=false)#check_length=false
sol_ns = solve(prob_ns,NewtonRaphson())

@variables z[1:2]
@parameters p[1:2]
ceq = [0 ~ z[1],0 ~ Symbolics.scalarize(cl[1])]
ns = NonlinearSystem(ceq,z,p)
z00=[2.2,2.3]
p00=[2.2,2.3]
prob_ns = NonlinearProblem(ns,z00,p00)
sol_ns = solve(prob_ns,NewtonRaphson())




nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
jac_func = generate_jacobian(ns)
f = @eval eval(nlsys_func)
Symbolics.jacobian(equations(ns),states(ns))
calculate_jacobian(ns)
f = eval(generate_function(ns, [x,y,z], [σ,ρ,β])[2])
cJ=generate_gradient(sys; parallel=Symbolics.MultithreadedForm())
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
jac_func = generate_jacobian(ns)
f = @eval eval(nlsys_func)
calculate_hessian

generate_hessian
hessian_sparsity
    jacobian_sparsity
    Joop, Jiip = eval.(generate_gradient(sys))
    
    multithreadedf = eval(ModelingToolkit.build_function(du,u,fillzeros=true,
                      parallel=ModelingToolkit.MultithreadedForm())[2])

#Joop(vars,params)
Joop(map(x->x[2],u0),map(x->x[2],p))
calculate_gradient(sys)

# define derivatives
f(x) = 2x + x^2
@register f(x)
function ModelingToolkit.derivative(::typeof(f), args::NTuple{1,Any}, ::Val{1})
    2 + 2args[1]
end
expand_derivatives(Dx(f(x)))

beta_ccdf(a,b,x) = zero(x) < x < one(x) ? one(x) - IncBetaDer.beta_inc_grad(a, b, x)[1] : zero(x)
function c_chance(z,p)
    c, α, β = selc*z, selα*z, selβ*z
    w, delta = p[1:N], p[2N+1:3N]
    zeroz = zero(z[1])+0.0001
    onez = one(z[1])-0.0001
    return beta_ccdf.(α,β,max.(min.(w./c,onez),zeroz) ) .- delta
end

## JuMP
using JuMP, Ipopt
m = Model(with_optimizer(Ipopt.Optimizer,
                            max_iter=100000,
                            print_level=5,
                            mu_strategy="adaptive",
                            derivative_test="second-order",
                            hessian_approximation="exact"
                        )
         )

# @variable(m, low[i] <= z[i=1:M] <= upp[i], start = z0[i])
# @objective(m, Max,sum(z[i]^2 for i=1:M))
# JuMP.optimize!(m)

## constraints and objective
# hᵀz=q
qJ, hᵀJ = p0[4N+2:6N+1], reshape(p0[6N+2:6N+1+2N*M],2N,M)
@constraint(m, hᵀJ*z .== qJ )

# zᵀHz = 0
HJ  = reshape(p0[6N+2+2N*M:6N+1+2N*M+M^2],M,M)
@constraint(m,dot(z,HJ*z) == 0 )

for i=1:N
    @eval ($(Symbol("c_chanceJ$i")))(z...) = c_chance(collect(z),p0)[$i]
    @eval JuMP.register(m, :($(Symbol("c_chanceJ$i"))) , M,  ($(Symbol("c_chanceJ$i"))), autodiff=true)
    @eval @NLconstraint(m,($(Symbol("c_chanceJ$i")))(z...) == 0)
end

objJ(z...)=obj_cont_int(collect(z),p0)
JuMP.register(m, :objJ, M, objJ, autodiff=true)
@NLobjective(m, Min , objJ(z...)  )

unset_silent(m)
JuMP.optimize!(m)

st0 = termination_status(m)
obj0 = objective_value(m)
@show st0, obj0

Asol0 = reshape(selA*JuMP.value.(z),N,N)
bsol0 = selb*JuMP.value.(z)
csol0 = selc*JuMP.value.(z)
αsol0 = selα*JuMP.value.(z)
βsol0 = selβ*JuMP.value.(z)
psol0 = reshape(selp*JuMP.value.(z),N,mini_batch)

zsol0 = [Asol0...;bsol0...;csol0...;psol0...]
Zsol0 = zsol0*zsol0'
rank(Zsol0)
eigvals(Zsol0)


cdf_pkg = [1-cdf(Beta(αsol0[i],βsol0[i]),w[i]/csol0[i]) for i=1:N]
cdf_opt = [dist_cdf([αsol0[i],βsol0[i],w[i]/csol0[i] ]...) for i=1:N]
tol = 1e-5
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol0.<=1)
    @test all(0 .<=bsol0.<=data.p_bar)
    @test all(0 .<=csol0.<=data.assets)
    @test all(0 .<=psol0.<=data.p_bar)

    @test norm(min.(p_bar,max.((1+g0).*(Asol0'*psol0[:,1] .+ csol0 .- x0[:,1].*csol0 ) .- g0.*p_bar,0))-psol0[:,1])<tol


    @test norm(cdf_pkg-cdf_opt)<0.01
    @test norm(delta-cdf_opt)<0.01

end


## NLP Models
using NLPModels, ADNLPModels, NLPModelsIpopt, ReverseDiff, NLPModelsModifiers
nc = length(_c(z0,p0) )
nlp = ADNLPModel(z->obj_cont_int(z,p0), z0, low, upp, z->_c(z,p0), zeros(nc), zeros(nc); lin=1:2N) #  adbackend = ADNLPModels.ReverseDiffAD()  adbackend = ADNLPModels.ZygoteAD() ; lin=1:2N)
obj(nlp,nlp.meta.x0)
gx = grad(nlp, nlp.meta.x0)
nlpLBFGS=NLPModelsModifiers.LBFGSModel(nlp)

function hess(nlp::AbstractNLPModel, x::AbstractVector; obj_weight::Real = one(eltype(x)))
  @lencheck nlp.meta.nvar x
  rows, cols = hess_structure(nlp)
  vals = hess_coord(nlp, x, obj_weight = obj_weight)
  Symmetric(sparse(FiniteDiff.finite_difference_hessian(x->obj(nlp,x),x)), :L)
end

hess(nlpLBFGS,z0)


Hx = hess(nlp, nlp.meta.x0)

obj(nlp, nlp.meta.x0) # evaluate the objective value at x0
grad!(nlp, nlp.meta.x0, g) # evaluate the objective gradient at x0
cons!(nlp, nlp.meta.x0, c) # evaluate the vector of constraints at x0


jac_structure!(nlp, rows, cols): fill rows and cols with the spartity structure of the Jacobian, if the problem is constrained;
jac_coord!(nlp, nlp.meta.x0, vals): fill vals with the Jacobian values corresponding to the sparsity structure returned by jac_structure!();
hess_structure!(nlp, rows, cols): fill rows and cols with the spartity structure of the lower triangle of the Hessian of the Lagrangian;
hess_coord!(nlp, x, y, vals; obj_weight=1.0): fill vals with the values of the Hessian of the Lagrangian corresponding to the sparsity structure returned by hess_structure!(), where obj_weight is the weight assigned to the objective, and y is the vector of multipliers.

stats = ipopt(nlpLBFGS)
print(stats)
problem = createProblem(n, nlp.meta.lvar, nlp.meta.uvar,
m, nlp.meta.lcon, nlp.meta.ucon,
nlp.meta.nnzj, nlp.meta.nnzh,
eval_f, eval_g, eval_grad_f,
eval_jac_g);#, eval_h)



using DCISolver, ADNLPModels, Logging


grad(nlp,z0) == ForwardDiff.gradient(z->obj_cont_int(z,p0),z0)
dci(nlp, nlp.meta.x0)

stats = with_logger(NullLogger()) do
  dci(nlp, nlp.meta.x0)
end

println(stats)




Pkg.add(url="https://github.com/JuliaSmoothOptimizers/NCL.jl", rev="master")
using NLPModelsIpopt, NCL
nlp = ADNLPModel(z->obj_cont_int(z,p0), z0, low, upp, z->_c(z,p0), zeros(nc), ones(nc)) # adbackend = ADNLPModels.ZygoteAD() ; lin=1:2N)

ncl_nlin_res = NCLModel(nlp)
NCLSolve(ncl_nlin_res)








# end #module



using PolyChaos
N=5
α=ones(N);β=2.0*ones(N);
mop = MultiOrthoPoly([ Beta01OrthoPoly(5, α[i],β[i]) for i in 1:N])







import Pkg; Pkg.add("SparseGrids");Pkg.add("FastGaussQuadrature");Pkg.add("PolyChaos")
using SparseGrids, FastGaussQuadrature, PolyChaos
#If f is a function that returns nodes, weights = f(n), for any (integer) order n, then the function sparsegrid computes the sparse extension to D dimensions of order O
nodes, weigths = sparsegrid(D, O, f)

order = [3, 5, 6, 4]
op1 = GaussOrthoPoly(order[1], addQuadrature=true);
op2 = Uniform01OrthoPoly(order[2], addQuadrature=true);
op3 = Beta01OrthoPoly(order[3], 2, 1.2, addQuadrature=true);
ops = [op1, op2, op3];
integrate(c->integrate(b->integrate((a)->a*b*c, op1),op2),op3)


order, Nrec, N, tol = 2, 3, 3, 1e-4
mop = MultiOrthoPoly([ Beta01OrthoPoly(order, α[i],β[i];  Nrec=Nrec) for i in 1:N], order)
α =rand(3) .+ 1.0; β=rand(3).+2.0
op1=Beta01OrthoPoly(order, α[1],β[1];  Nrec=Nrec)
op2=Beta01OrthoPoly(order, α[2],β[2];  Nrec=Nrec)
op2=Beta01OrthoPoly(order, α[3],β[3];  Nrec=Nrec)
function f()
    integrate(c->integrate(b->integrate((a)->a, op1),op2),op3)
end
α[1]/(α[1]+β[1])
@btime $f()

show(mop)
calculateAffinePCE(mop)
Quad(10,mop.measure)
a,b = coeffs(mop)
mean(coeffs(mop),mop)
calculateAffinePCE(1,mop)

# all-purpose constructor (last resort!)
function Quad(N::Int, w::Function, dom::Tuple{<:Real,<:Real}; quadrature::Function=clenshaw_curtis)
    N <= 0 && throw(DomainError(N, "number of quadrature points has to be positive"))
    nodes, weights = quadgp(w, dom[1], dom[2], N; quadrature=quadrature)
    Quad("quadgp", N, nodes, weights)
end

function Quad(N::Int, measure::AbstractMeasure; quadrature::Function=clenshaw_curtis)
    typeof(measure) != Measure && @warn "For measures of type $(typeof(measure)) the quadrature rule should be based on the recurrence coefficients."
    Quad(N, measure.w, (measure.dom[1], measure.dom[2]); quadrature=quadrature)
end







Tensor(2,mop)

PolyChaos.integrate(x-> x[1]+x[2], mop)

nodes, weights = sparsegrid(3,20,n->nw(n,1.0,2.0))
# Evaluate integrand
F = map(x -> x, nodes)
Fn = map(sum, F)

# Evaluate quadrature sum
Q = dot(Fn, weights)



using PolyChaos, Test, LinearAlgebra
import QuadGK: quadgk

function integration_test(op::AbstractOrthoPoly, dim::Int,name::String;tol::Float64=1e-6,
                          dom::Tuple{<:Real,<:Real}=op.measure.dom)
    N, w, ind = op.deg, op.measure.w, zeros(Int64,dim)
    @testset "$name" begin
        for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
            ind[:] .= ind_[:]
            s1 = computeSP(ind, op)
            f(t) = prod(evaluate(ind[i], t, op) for i in 1:dim)*w(t)
            s2 = quadgk(f, dom[1], dom[2])[1]
            if abs(s1)<=1e-3
                @test isapprox(s1, s2; atol=tol)
            else
                @test isapprox(abs(s1/s2), 1; atol=tol)
            end
        end
    end
end

N, Nrec, dim, tol = 3, 1000, 3, 1e-4
α, β = 1.32, 3.23
op = Beta01OrthoPoly(N, α, β; Nrec=Nrec)
set_beta = integration_test(op, dim, "beta01"; tol=tol)
N, w, ind = op.deg, op.measure.w, zeros(Int64,dim)
ind_ = (1,2,3)
dom = (0,1)

for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
    ind[:] .= ind_[:]
    s1 = computeSP(ind, op)
    f(t) = prod(evaluate(ind[i], t, op) for i in 1:dim)*w(t)
    s2 = quadgk(f, dom[1], dom[2])[1]
    if abs(s1)<=1e-3
        @test isapprox(s1, s2; atol=tol)
    else
        @test isapprox(abs(s1/s2), 1; atol=tol)
    end
end

for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
   display(ind_[:] )
end



using PolyChaos
dd = [3]
op1 = JacobiOrthoPoly(dd[1],1.1,2.3)
op2 = JacobiOrthoPoly(dd[1],1.1,2.3)
show(op1)
ops = [op1,op2]
mop = MultiOrthoPoly(ops, dd[1])
mop.ind

using Revise
using SmolyakApprox, FastGaussQuadrature, LinearAlgebra
d  = 1  # Set the number of dimensions
mu = 3  # Set the level of approximation
g1, m1 = smolyak_grid(chebyshev_gauss_lobatto,d,mu)  # Construct the Smolyak grid and the multi index
g2, m2 =  smolyak_grid(clenshaw_curtis_equidistant,d,mu)
gj(n::Integer) = gaussjacobi(n,1.1,2.3)
g3, m3 = smolyak_grid(gj,d,mu)

function test(grid)
    #y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,3].+2).^0.2.*(grid[:,4].+2).^0.8.*(grid[:,5].+7).^0.1
    #y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,1].+2).^0.2.*(grid[:,2].+2).^0.8.*(grid[:,1].+7).^0.1
    y_value = (grid[:,1].+1).^0.1
    return y_value
end
y1 = test(g1)
y2 = test(g2)
y3 = test(g3)

w1 = smolyak_weights(y1,g1,m1)
w2 = smolyak_weights(y2,g2,m2)
w3 = smolyak_weights(y3,g3,m3)


im1 = smolyak_inverse_interpolation_matrix(g1,m1)  # Compute the interpolation matrix
im2 = smolyak_inverse_interpolation_matrix(g2,m2)  # Compute the interpolation matrix
im3 = smolyak_inverse_interpolation_matrix(g3,m3)  # Compute the interpolation matrix

w1 = smolyak_weights(y1,im1)             # Compute the Smolyak weights
w2 = smolyak_weights(y2,im2)             # Compute the Smolyak weights
w3 = smolyak_weights(y3,im3)             # Compute the Smolyak weights

#point = [0.75, 0.45, 0.82, -0.15, -0.95]
point = [0.75]
y_hatt1 = smolyak_evaluate(w1,point,m1)  # Evaluate the approximated function
y_hatt2 = smolyak_evaluate(w2,point,m2)  # Evaluate the approximated function
y_hatt3 = smolyak_evaluate(w3,point,m3)  # Evaluate the approximated function

  multi_index        = SmolyakApprox.generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = SmolyakApprox.m_i(unique_multi_index)

T= Float64
  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_nodes_gj   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights_gj = Array{Array{T,1},1}(undef,length(unique_node_number))
  for i = 1:length(unique_node_number)
    base_nodes_gj[i], base_weights_gj[i] = gj(unique_node_number[i])
    base_nodes[i], base_weights[i] = chebyshev_gauss_lobatto(unique_node_number[i])
  end

  unique_base_nodes_gj = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes_gj[1] = base_nodes_gj[1]
  unique_base_nodes = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes[1] = base_nodes[1]
  for i = 2:length(unique_base_nodes)
    unique_base_nodes[i] = setdiff(base_nodes[i],base_nodes[i-1])
    unique_base_nodes_gj[i] = setdiff(base_nodes_gj[i],base_nodes_gj[i-1])
  end


  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  #nodes_gj = Array{T,2}(undef,determine_grid_size(multi_index))
  nodes_gj = Array{T,2}(undef,(sum(map(x->length(x),unique_base_nodes_gj)), d))
  l = 1
  l_gj = 1
  @inbounds for j = 1:size(multi_index,1)
    new_nodes = unique_base_nodes[multi_index[j,1]]  # Here new_nodes is a 1d array
    new_nodes_gj = unique_base_nodes_gj[multi_index[j,1]]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[j,i]])  # Here new_nodes becomes a 2d array
      new_nodes_gj = combine_nodes(new_nodes_gj,unique_base_nodes_gj[multi_index[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    m_gj = size(new_nodes_gj,1)
    nodes[l:l+m-1,:] = new_nodes
    nodes_gj[l_gj:l_gj+m_gj-1,:] = new_nodes_gj
    l += m
    l_gj += m_gj
  end






function determine_grid_size(mi)

  temp = similar(mi)

  for i = 1:size(mi,1)
    for j = 1:size(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      elseif mi[i,j] == 2
        temp[i,j] = 2^(mi[i,j]-1)
      else
        temp[i,j] = 2^(mi[i,j]-1)+1 - (2^(mi[i,j]-2)+1)
      end
    end
  end

  s = 0
  for i = 1:size(mi,1)
    t = 1
    for j = 1:size(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end
#= A second way of computing the weights and evaluating the approximated function that
  computes the interpolation matrix just once. =#

  interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix
  w = smolyak_weights(y,interp_mat)             # Compute the Smolyak weights
  y_hatt = smolyak_evaluate(w,point,multi_ind)  # Evaluate the approximated function

  # Piecewise linear
 y_pl = test(g2)
  w_pl = smolyak_pl_weights(y_pl,g2,m2)
  y_pl_hat = smolyak_pl_evaluate(w_pl,point,g2,m2)

  # Evaluate the exact function at point

  y_actual = test(point')
y_hatt


function nw(n::S,domain = [1.0,-1.0])  where {S<:Integer}
    a=1.0
    b=2.0
    opq = Beta01OrthoPoly(n, a, b, addQuadrature = true)
    opq = Beta01OrthoPoly(n, a, b, addQuadrature = true)
    jacobi
    return (opq.quad.nodes,opq.quad.weights)
end

clenshaw_curtis_equidistant(5)
nw(5)

y = grid[:,1]
w = smolyak_weights(y,interp_mat)

interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix


function chebyshev_gauss_lobatto(n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # These nodes a just the Chebyshev extrema by a different name.

  # Construct the nodes on the [-1.0,1.0] interval

  if n <= 0
    error("The number of nodes must be positive.")
  end

  if n == 1
    nodes   = [0.0]
    weights = [1.0*pi]
  else
    nodes    = zeros(n)
    nodes[1] = -1.0
    nodes[n] = 1.0

    weights    = zeros(n)
    weights[1] = pi/(2*(n-1))
    weights[n] = pi/(2*(n-1))

    for i = 2:div(n,2)
      x = cos(pi*(i-1)/(n-1))
      nodes[i]       = -x
      nodes[end-i+1] = x
      weights[i]       = pi/(n-1)
      weights[end-i+1] = pi/(n-1)
    end

    if isodd(n)
      nodes[round(Int,(n+1)/2)]   = 0.0
      weights[round(Int,(n+1)/2)] = pi/(n-1)
    end
  end

  # Scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes, weights

end







using Surrogates, PolyChaos
using Plots
n = 20
lower_bound = 1.0
upper_bound = 6.0
x = sample(n,lower_bound,upper_bound,LowDiscrepancySample(2))
f = x -> log(x)*x + sin(x)
y = f.(x)
scatter(x, y, label="Sampled points", xlims=(lower_bound, upper_bound), legend=:top)
plot!(f, label="True function", xlims=(lower_bound, upper_bound), legend=:top)

poly1 = PolynomialChaosSurrogate(x,y,lower_bound,upper_bound)
poly2 = PolynomialChaosSurrogate(x,y,lower_bound,upper_bound, op = GaussOrthoPoly(5))
plot(x, y, seriestype=:scatter, label="Sampled points", xlims=(lower_bound, upper_bound), legend=:top)
plot!(f, label="True function",  xlims=(lower_bound, upper_bound), legend=:top)
plot!(poly1, label="First polynomial",  xlims=(lower_bound, upper_bound), legend=:top)
plot!(poly2, label="Second polynomial",  xlims=(lower_bound, upper_bound), legend=:top)




## test quadrature
using Quadrature, ForwardDiff, FiniteDiff, Zygote, Cuba
f(x,p) = sum(sin.(x .* p))
lb = ones(2)
ub = 3ones(2)
p = [1.5,2.0]

function testf(p)
    prob = QuadratureProblem(f,lb,ub,p)
    sin(solve(prob,CubaCuhre(),reltol=1e-6,abstol=1e-6)[1])
end
dp1 = Zygote.gradient(testf,p)
dp2 = FiniteDiff.finite_difference_gradient(testf,p)
dp3 = ForwardDiff.gradient(testf,p)
dp1[1] ≈ dp2 ≈ dp3



    
    
using ModelingToolkit, SparseArrays, Test, GalacticOptim, Optim

@variables x y
@parameters a b
f(p1,p2)= p1+a*p2
f(x,y)
    
loss = (a - x)^2 + b * (y - x^2)^2
sys1 = OptimizationSystem(loss,[x,y],[a,b],name=:sys1)
sys2 = OptimizationSystem(loss,[x,y],[a,b],name=:sys2)

@variables z
@parameters β
loss2 = sys1.x - sys2.y + z*β
combinedsys = OptimizationSystem(loss2,[z],[β],systems=[sys1,sys2],name=:combinedsys)

equations(combinedsys)
states(combinedsys)
parameters(combinedsys)

calculate_gradient(combinedsys)
calculate_hessian(combinedsys)
generate_function(combinedsys)
generate_gradient(combinedsys)
generate_hessian(combinedsys)
ModelingToolkit.hessian_sparsity(combinedsys)
    
using ModelingToolkit  ,Symbolics
@variables x[1:2] y[1:3] A[1:3,1:3]
@parameters σ[1:2]
eqs = [0 ~ σ[1]*x[1]*y[1]+σ[2]*x[2]+y[2]*y[3]]
ns = NonlinearSystem(eqs, vcat(x...,y...,A...),σ)
calculate_jacobian(ns)
generate_jacobian(ns)
generate_hessian(ns)
ModelingToolkit.jacobian_sparsity(ns)
structural_simplify(ns)
ns
    
    
