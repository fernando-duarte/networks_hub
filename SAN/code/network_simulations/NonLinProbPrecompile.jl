
# module NonLinProbPrecompile
# include("NetDefs.jl"); using .NetDefs
# using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
# export f_noeval_good
# function system(; kwargs...)
#     # Define some variables
#     @variables z[1:M]
#     @parameters p[1:P]
#     zcat = vcat(z...);
#     pcat = vcat(p...);
#     # Define a system of nonlinear equations
#     ceq = vcat(
#         0 .~ c_lin(zcat,pcat),
#         0  ~ c_quad(zcat,pcat),
#         #0 .~ c_p(zcat,pcat,x0),
#         0 .~ c_chance(zcat,pcat)
#       )
#     ns = NonlinearSystem(ceq,z,p)
#     return generate_function(ns,z,p)
# end
# # Setting eval_expression=false and eval_module=[this module] will ensure
# # the RGFs are put into our own cache, initialised below.   
# RuntimeGeneratedFunctions.init(@__MODULE__)
# const f_noeval_good = system(; eval_expression=false, eval_module=@__MODULE__)
# end #module

# module NonLinProbPrecompileObj
#     #include("NetworkType.jl");using .NetworkType
#     include("NetDefs.jl"); using .NetDefs
#     using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
#     export f_noeval_obj
#     function system(; kwargs...)
#         # Define some variables
#         @variables z[1:M]
#         zcat = vcat(z...);
#         # Define a system of nonlinear equations
#         p0A = Array(p0)
#         loss = obj(zcat,p0A)
#         opt = OptimizationSystem(loss,z,[])
#         return generate_function(opt,z,[])
#     end
#     # Setting eval_expression=false and eval_module=[this module] will ensure
#     # the RGFs are put into our own cache, initialised below.   
#     RuntimeGeneratedFunctions.init(@__MODULE__)
#     const f_noeval_obj = system(; eval_expression=false, eval_module=@__MODULE__)
# end #module


module NonLinProbPrecompileNum
include("NetDefs.jl"); using .NetDefs
using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
export f_noeval_num
function system(; kwargs...)
    # Define some variables
    @variables z[1:M]
    zcat = vcat(z...);
    # Define a system of nonlinear equations
    p0A = Array(p0)
    ceq = vcat(
        0 .~ c_lin(zcat,p0A),
        0  ~ c_quad(zcat,p0A),
        #0 .~ c_p(zcat,p0A,x0),
        0 .~ c_chance(zcat,p0A)
      )
    ns = NonlinearSystem(ceq,z,[])
    return generate_function(ns,z,[])
end
# Setting eval_expression=false and eval_module=[this module] will ensure
# the RGFs are put into our own cache, initialised below.   
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_num = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module NonLinProbPrecompileContractionQuadUnif
include("NetDefs.jl"); using .NetDefs
using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
export f_noeval_contraction_quaduniform
function system(; kwargs...)
    # Define some variables
    @variables z[1:M]
    zcat = vcat(z...);
    p0A = Array(p0)
    cont = obj_contraction_quaduniform(zcat, p0A, x0,weights0)  
    opt = OptimizationSystem(cont,z,[])
    return generate_function(opt,z,[])
end
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_contraction_quaduniform = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module NonLinProbPrecompileSumpContraction
include("NetDefs.jl"); using .NetDefs
using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
export f_noeval_sump_contraction
function system(; kwargs...)
    # Define some variables
    @variables z[1:M]
    zcat = vcat(z...);
    p0A = Array(p0)
    cont = sump_contraction(zcat, p0A, x0)
    opt = OptimizationSystem(cont,z,[])
    return generate_function(opt,z,[])
end
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_sump_contraction = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module PrecompileExtpow
include("NetDefs.jl"); using .NetDefs
using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
export f_noeval_extpow_prod_pdf
function system(; kwargs...)
    # Define some variables
    @variables z[1:M]
    zcat = vcat(z...);
    cont = extpow_prod_pdf(x0,zcat) 
    opt = OptimizationSystem(cont,z,[])
    return generate_function(opt,z,[])
end
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_extpow_prod_pdf = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module NonLinProbPrecompileObj
include("NetDefs.jl"); using .NetDefs
using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
export f_noeval_obj
function system(; kwargs...)
    # Define some variables
    @variables z[1:M]
    zcat = vcat(z...);
    p0A = Array(p0)
    cont = obj_contraction(zcat, p0A, x0)
    opt = OptimizationSystem(cont,z,[])
    return generate_function(opt,z,[])
end
RuntimeGeneratedFunctions.init(@__MODULE__)
const f_noeval_obj = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

