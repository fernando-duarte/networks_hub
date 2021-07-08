module NonLinProbPrecompile
    #include("NetworkType.jl");using .NetworkType
    include("NetDefs.jl"); using .NetDefs
    using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
    export f_noeval_good
    function system(; kwargs...)
        # Define some variables
        @variables z[1:M]
        @parameters p[1:P]
        zcat = vcat(z...);
        pcat = vcat(p...);
        # Define a system of nonlinear equations
        ceq = vcat(
            0 .~ c_lin(zcat,pcat),
            0  ~ c_quad(zcat,pcat),
            0 .~ c_p(zcat,pcat,x0)
          )
        ns = NonlinearSystem(ceq,z,p)
        return generate_function(ns,z,p)
    end
    # Setting eval_expression=false and eval_module=[this module] will ensure
    # the RGFs are put into our own cache, initialised below.   
    RuntimeGeneratedFunctions.init(@__MODULE__)
    const f_noeval_good = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module NonLinProbPrecompileNum
    #include("NetworkType.jl");using .NetworkType
    include("NetDefs.jl"); using .NetDefs
    using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
    export f_noeval_num
    function system(; kwargs...)
        # Define some variables
        @variables z[1:M]
        zcat = vcat(z...);
        # Define a system of nonlinear equations
        ceq = vcat(
            0 .~ c_lin(zcat,p_dense),
            0  ~ c_quad(zcat,p_dense),
            0 .~ c_p(zcat,p_dense,x0)
          )
        ns = NonlinearSystem(ceq,z,[])
        return generate_function(ns,z,[])
    end
    # Setting eval_expression=false and eval_module=[this module] will ensure
    # the RGFs are put into our own cache, initialised below.   
    RuntimeGeneratedFunctions.init(@__MODULE__)
    const f_noeval_num = system(; eval_expression=false, eval_module=@__MODULE__)
end #module

module NonLinProbPrecompileObj
    #include("NetworkType.jl");using .NetworkType
    include("NetDefs.jl"); using .NetDefs
    using ModelingToolkit, LinearAlgebra, RuntimeGeneratedFunctions
    export f_noeval_obj
    function system(; kwargs...)
        # Define some variables
        @variables z[1:M]
        zcat = vcat(z...);
        # Define a system of nonlinear equations
        loss = obj(zcat,p_dense)
        opt = OptimizationSystem(loss,z,[])
        return generate_function(opt,z,[])
    end
    # Setting eval_expression=false and eval_module=[this module] will ensure
    # the RGFs are put into our own cache, initialised below.   
    RuntimeGeneratedFunctions.init(@__MODULE__)
    const f_noeval_obj = system(; eval_expression=false, eval_module=@__MODULE__)
end #module
