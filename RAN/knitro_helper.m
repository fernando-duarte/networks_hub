
function [xstar,fval,exitflag,output,lambda,grad, t1] = knitro_helper(x0,Aineq,bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, p_bar, c, b, w, g ,delta, N, num_sim, use_parallel, index_nan_b, index_nan_c, uncorr_shocks)
    tester = 0;
    stopper = 0;
    j = 1;
    while stopper == 0
        try
        tester = tester + 1;
        tic
        obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),b(~index_nan_b),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c, -1, 0, uncorr_shocks, 'uniform'); 
        [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),x0,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, file_name);%,'knitro.opt'
        t1 = toc;
        stopper = 1;s
        catch err
            disp(tester)
        end
    end
end