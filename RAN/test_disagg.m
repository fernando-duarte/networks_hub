%%

% cd('/home/rcecdj02/Fernando/Networks_Simulation_Local/Matlab')
clear
clc
clear simulation_ratio

%ORDER IN VECTOR: Astar - b (outside liabilities) - c (outside assets)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Parameters for Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configurations for run -- Plotting
%Tests of numeric-vs-analytic gradient
gradcheck = false;
%Whether to make pdf plots of post-sim networks
plot_network = true;
%Whether to make pdf plot of objective value and NVI
plot_nvi = false;
%Whether to plot series of sims with diff numbers of uncorrelated,
%underlying shocks
plot_uncorrs = false;

plot_obj_fun_test = false;

run_analytics_tests = false;
run_nvi_tightness = false;

benchmark = true;
%Configurations for run -- Simulations to run
%Test diff numbers of uncorrelated underlying shocks
uncorr_shock_test = false;
%Run multistart (options in knitro_me_multistart.opt)
multistart = false;
%Run knitro again, at endpoint of first sim
rerun = false;
%Test objective function value at endpoint (run with 10x num_sim)
obj_fun_test = false;

use_parallel = 0; 
pool_size = 8;
%Dissaggregate aggregate nodes
diss_agg = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Network characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Range of quarters to run. 192 = 2008Q1
    %194 is a generally interesting quarter
q_start = 194;
q_end = 194;

%Number of simulations to calculate objective function
num_sim = 10;

%Num of nodes in network (will take top N by assets)
N = 3;

%Default number of uncorrelated underlying shocks
uncorr_shocks = N;

%Bankruptcy costs
g = 0;

missing_beta_ineq = 1;
missing_beta_ineq_val = .90;

addpath(genpath('Toolboxes'))

% sheets ={'Broker Dealers', 'Insurance','Other','REITs','BHCs'};
% nsheets = length(sheets);

sort_by = 'assets';
spec_id = [datestr(date) '_numsim' num2str(num_sim) '_top' num2str(N) 'by' sort_by '_qstart' num2str(q_start) '_qend' num2str(q_end) '_betamax' num2str(100*missing_beta_ineq_val) 'perc'];

%Diary will be saved here
diaryfile = ['./Diaries/sim_results_' spec_id];
%Output data will be saved here
to_save_file = ['./Simulation_Results/sim_results_' spec_id];
diary(diaryfile)
diary on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Building structs for saved results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Benchmark run
results_save_worst = struct('qt', {}, 'options', {}, 'exit', {}, 'times', {}, ...
    'output', {}, 'options_file', {}, 'output_value', {}, 'obj_value', {}, ...
    'w', {}, 'a', {}, 'names', {}, 'p_bar', {}, 'c', {}, 'delta', {}, 'nan_c', {}, ...
    'nan_b', {}, 'nvi', {}, 'betas', {}, 'sectors', {}, 'tkr', {}, 'nvi_benchmark', {}, ...
    'uncorr_shocks', {}, 'Astar', {});

%Tests of uncorrelated shocks
results_save_worst_uncorr_iter = struct('qt', {}, 'options', {}, 'exit', {}, 'times', {}, ...
    'output', {}, 'options_file', {}, 'output_vresults_save_worst_uncorr_iteralue', {}, 'obj_value', {}, ...
    'w', {}, 'a', {}, 'names', {}, 'p_bar', {}, 'c', {}, 'delta', {}, 'nan_c', {}, 'nan_b', {}, 'nvi', {}, 'betas', {}, 'sectors', {}, 'tkr', {}, 'nvi_benchmark', {}, 'uncorr_shocks', {});

%Second runthrough
results_save_worst2 = struct('qt', {}, 'options', {}, 'exit', {}, 'times', {}, ...
    'output', {}, 'options_file', {}, 'output_value', {}, 'obj_value', {}, ...
    'w', {}, 'a', {}, 'names', {}, 'p_bar', {}, 'c', {}, 'delta', {}, 'nan_c', {}, 'nan_b', {}, 'nvi', {}, 'betas', {}, 'sectors', {}, 'tkr', {}, 'nvi_benchmark', {}, 'uncorr_shocks', {});

%Test of multistart
results_save_worst_multi = struct('qt', {}, 'options', {}, 'exit', {}, 'times', {}, ...
    'output', {}, 'options_file', {}, 'output_value', {}, 'obj_value', {}, ...
    'w', {}, 'a', {}, 'names', {}, 'p_bar', {}, 'c', {}, 'delta', {}, 'nan_c', {}, ...
    'nan_b', {}, 'nvi', {}, 'betas', {}, 'sectors', {}, 'tkr', {}, 'nvi_benchmark', {}, ...
    'uncorr_shocks', {}, 'Astar', {});
%Best-case network
results_save_best = struct('qt', {}, 'options', {}, 'exit', {}, 'times', {}, ...
    'output', {}, 'options_file', {}, 'output_value', {}, 'obj_value', {}, ...
    'w', {}, 'a', {}, 'names', {}, 'p_bar', {}, 'c', {}, 'delta', {}, 'nan_c', {}, 'nan_b', {}, 'nvi', {}, 'betas', {}, 'sectors', {}, 'tkr', {}, 'nvi_benchmark', {}, 'uncorr_shocks', {});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Starting pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_parallel
    p = gcp('nocreate'); 
    if isempty(p)
         % if pool not open, create one
         mypool=parpool(pool_size,'AttachedFiles',{'ui.m','systemeq.m','sums.m'});
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Running chosen simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for q=q_start:q_end
    disp('-----------------STARTING----------------------')
    disp(q)

    load('network_data_all')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Selecting, cleaning data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Converting from cell arrays
    scale = 1e6;
    p_bar = cell2mat(p_bar)/scale;
    c = cell2mat(c)/scale;
    b = cell2mat(b)/scale;
    w = cell2mat(w)/scale;
    delta = cell2mat(delta);
    a = cell2mat(a)/scale;
    d = cell2mat(d)/scale;
    f = cell2mat(f)/scale;
    qt = cell2mat(qt);
    beta = cell2mat(beta);
    type = cell2mat(type);
    names = vertcat(names{:});
    sectors = vertcat(sectors{:});
    tkr = vertcat(tkr{:});
    nvi_benchmark = cell2mat(nvi);

    %Isolating this loop's quarter
    right_qt = qt == q;
    p_bar = p_bar(right_qt);
    c = c(right_qt);
    b = b(right_qt);
    w = w(right_qt);
    delta = delta(right_qt);
    a = a(right_qt);
    d = d(right_qt);
    f = f(right_qt);
    type = type(right_qt);
    names = names(right_qt);
    tkr = tkr(right_qt);
    beta = beta(right_qt);
    sectors = sectors(right_qt);
    nvi_benchmark = nvi_benchmark(right_qt);
    
    %Tossing nodes without total assets data
    keep = ~isnan(a);
    a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep); beta=beta(keep); %index_top10 =index_top10(keep); 
    names = names(keep); sectors = sectors(keep); nvi_benchmark = nvi_benchmark(keep);
    
    %Tossing nodes with negative net worth
    keep = w>0;
    a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep); beta = beta(keep); %index_top10 =index_top10(keep); 
    names = names(keep); sectors = sectors(keep); nvi_benchmark = nvi_benchmark(keep); tkr = tkr(keep);
   
    %Only keeping nodes that could conceivably default on their own (not a
    %big deal)
    keep = c > w |  isnan(c);
    a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep); beta = beta(keep); %index_top10 =index_top10(keep); 
    names = names(keep); sectors = sectors(keep); nvi_benchmark = nvi_benchmark(keep); tkr = tkr(keep);

    %Sort by assets, take top N for simulation
    if strcmp(sort_by, 'assets')
        [sorted,ix] = sortrows(a,-1);
    elseif strcmp(sort_by, 'delta')
        [sorted,ix] = sortrows(delta,-1);        
    end
    
    a = a(ix); c = c(ix);b= b(ix);w = w(ix);delta= delta(ix);d = d(ix);f = f(ix);type = type(ix);p_bar =p_bar(ix); beta=beta(ix); %index_top10 =index_top10(ix); 
    names = names(ix); sectors = sectors(ix); nvi_benchmark = nvi_benchmark(ix); tkr = tkr(ix);
    a = a(1:N)';c = c(1:N)';b= b(1:N)';w = w(1:N)';delta= delta(1:N)';d = d(1:N)';f = f(1:N)';type = type(1:N)';p_bar =p_bar(1:N)'; beta=beta(1:N)'; %index_top10 =index_top10(1:N)'; 
    names = names(1:N)'; sectors = sectors(1:N)'; nvi_benchmark = nvi_benchmark(1:N)'; tkr = tkr(1:N)';

    %REMOVE THIS -- just for testing
    %order: oth - ins - JPM
        b = [1,2,1]; 
        c = [3,7,2];
       
%           b = [2,1,1]; 
%           c = [7,3,2];
       
%        d = [5, NaN, 5];
%        f = [5, NaN, 5];       
       
       
       d = [5, 5, 5];
       f = [5, 5, 5];
        
       p_bar = f+b;
       a = c+d;
       w= a-p_bar;
       
%      f = nan(size(f)); %this seems to be crashing the program
%      d = nan(size(d));  %this causes an error in plotting
%      beta=nan(size(beta));    
    
    %We are including both aggregate approx sector nodes, and (potentially)
    %individual nodes from those sectors. Deduct the assets of included nodes
    %from the included assets of the approx sector node. Preserve the original
    %ratio of inside/outside assets (which was determined from the Flow of
    %Funds)
    insur = find(strcmp(sectors,'Insurance'));
    other = find(strcmp(sectors,'Other'));
    reit  = find(strcmp(sectors,'REITs'));

    insur_agg = find(strcmp(names,'Insurance Aggregate'));
    reit_agg  = find(strcmp(names,'REIT Aggregate'));
    other_agg = find(strcmp(names,'Other Aggregate'));

    c_ratio_insur = c(insur_agg)/a(insur_agg);
    c_ratio_other = c(other_agg)/a(other_agg);
    c_ratio_reit  = c(reit_agg)/a(reit_agg);
    %uncomment when running the main code
%     a(insur_agg)     = a(insur_agg) - sum(a(insur));
%     p_bar(insur_agg) = p_bar(insur_agg) - sum(p_bar(insur));
%     w(insur_agg)     = a(insur_agg) - p_bar(insur_agg);
%     c(insur_agg)     = c_ratio_insur * a(insur_agg);
%     
%     a(other_agg)     = a(other_agg) - sum(a(other));
%     c(other_agg)     = c_ratio_other * a(other_agg);
%     p_bar(other_agg) = p_bar(other_agg) - sum(p_bar(other));
%     w(other_agg)     = a(other_agg) - p_bar(other_agg);
%     
%     a(reit_agg)      = a(reit_agg) - sum(a(reit));
%     c(reit_agg)      = c_ratio_reit * a(reit_agg);
%     p_bar(reit_agg)  = p_bar(reit_agg) - sum(p_bar(reit));
%     w(reit_agg)      = a(reit_agg) - p_bar(reit_agg);
    
    %Calculate NVI with simulation nodes (at appropriate gamma)
    beta_plus = max(beta);
    delta_a   = a .* delta;
    nvi = sum(delta_a)/((1-(1+g)*beta_plus)*sum(a));

    % find missing values for b and c
    index_nan_b = isnan(b);
    index_nan_c = isnan(c);
    index_nan   = index_nan_b | index_nan_c;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Building constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Constraints
    if missing_beta_ineq
    Aineq = sparse(repmat(eye(N),1,N)); bineq=a;
    %Aineq=repmat(eye(N),1,N); bineq = (f./p_bar)'; % sum of columns of A cannot exceed ratio of inside to total liabilities
    Aineq = [Aineq,sparse(N,sum(index_nan_b)),sparse(N,sum(index_nan_c))];
    bineq(index_nan_c) = missing_beta_ineq_val;
    else
       Aineq = []; bineq = []; 
    end
        
    % linear equality constraints for non-missing b and c
    Aeq1=sparse(kron(eye(N),p_bar)); beq1=d';  % for each j, sum of p_ij across i = inside assets of j
    Aeq2=sparse(repmat(diag(p_bar),1,N)); beq2 =f';  % for each i, sum of p_ij across j = inside liabilities of i
    Aeq3=sparse(1:N,sub2ind([N,N],1:N,1:N),1,N,N^2);beq3=sparse(N,1); % diagonal elements of A must be zero
    
    % linear inequality constraint adding missing b and c
    temp = diag(index_nan_c);
    temp(:,all(temp==0))=[];
    Aeq1 = [Aeq1,sparse(N,sum(index_nan_b)),temp]; beq1(index_nan_c)=a(index_nan_c);

    temp = diag(index_nan_b);
    temp(:,all(temp==0))=[];
    Aeq2 = [Aeq2,temp,sparse(N,sum(index_nan_c))]; beq2(index_nan_b)=p_bar(index_nan_b);

    Aeq3 = [Aeq3,sparse(N,sum(index_nan_b)),sparse(N,sum(index_nan_c))];

    Aeq = [Aeq1;Aeq2;Aeq3];
    beq = [beq1;beq2;beq3]; 
  
    %Lower, upper bounds
    lb_A=sparse(N^2,1);   
    ub_A = repmat((f./p_bar)',N,1); 
    ub_A(isnan(ub_A))=1;
    lb_b=sparse(sum(index_nan_b),1)+1e-5;   ub_b = ones(sum(index_nan_b),1).*p_bar(index_nan_b)'- 1e-5;
    lb_c=ones(sum(index_nan_c),1).*w(index_nan_c)'+1e-5;   ub_c = ones(sum(index_nan_c),1).*a(index_nan_c)'-1e-5;
    lb = [lb_A;lb_b;lb_c];
    ub = [ub_A;ub_b;ub_c];
    
    %replace above with dissagg constraints if specified
    if diss_agg
        disagg;
    end
    
    %Complementarity constraints
    nonlcon = @(x) sums(x,N,sum(index_nan_b)+sum(index_nan_c));

    upper_triang = ~tril(ones(N),0);
    lower_triang = ~triu(ones(N),0);
    [ix_row,ix_column]=find(upper_triang);
    ix_transpose = fliplr([ix_row ix_column]);
    ix =sortrows([sub2ind([N N],ix_transpose(:,1),ix_transpose(:,2)) find(upper_triang)],2);

    extendedFeatures.ccIndexList1 = ix(:,1)';
    extendedFeatures.ccIndexList2 = ix(:,2)';

    %initial point
    x0 = lb+(ub-lb)/2;
    if diss_agg
        x0=x0helper;
    end
    %x0 = [0;.08;.57;.25;0;0;0;.14;0;4;4]; 
    %x0 = [0;.04;.3;.1;0;0;0;.08;0;4;4]; 
    %x0 = [0.7090;0.8071;0.4846;0.7878;0.1892;0.4215;0.0634;0.5243;0.2812;2.4081;3.7838];
    iterations_to_run = 1;

    options = optimset('Display','iter-detailed','TolX',1e-4,'TolFun',1e-2, ...
        'TolCon', 1e-3,'GradConstr','on', 'GradObj','on');
    
    %use numerical gradient for disaggregated setup 
    if diss_agg
        options = optimset('Display','iter-detailed','TolX',1e-4,'TolFun',1e-2, ...
        'TolCon', 1e-3,'GradConstr','on', 'GradObj','off'); %this was edited
    end
    
    % Some scratch work from GlobalSearch testing
    % opts_globsearch = optimoptions(@fmincon, 'UseParallel', false, 'MaxFunctionEvaluations', 10000)
    % obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,1000,use_parallel,index_nan_b,index_nan_c, -1, 1); 
    % prob = createOptimProblem('fmincon', 'objective', obj, 'x0', x0, 'Aeq', Aeq, 'beq', beq, 'Aineq', Aineq, 'bineq', bineq, 'lb', lb, 'ub', ub, 'nonlcon', nonlcon, 'options', opts_globsearch);
    % g = GlobalSearch('OutputFcn', 'TestOutputFun', 'Display', 'iter', 'NumTrialPoints', 151, 'NumStageOnePoints', 150);
    % run(g, prob)

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Simulation runs
    %       Note: file_name is the knitro options file for each
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_name = char(strcat('knitro_benchmark.opt'));

    %Benchmark
    if benchmark
            [xstar,fval,exitflag,output,lambda,grad,t1] = knitro_helper(x0,Aineq,...
                bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, p_bar, c, b, ...
                w, g ,delta, N, num_sim, use_parallel, index_nan_b, index_nan_c, ...
                uncorr_shocks);
            Astar_ratio_worst = reshape(xstar(1:N^2, 1),N,N);
            if diss_agg
                disagg_cleanout
            end
            results_temp = struct('qt', {q}, 'options', {fileread(file_name)}, ...
                'exit', {exitflag}, 'times', {t1}, 'output', {xstar}, ...
                'options_file', {file_name}, 'output_value', {output}, 'obj_value', {fval}, ...
                'w', {w}, 'a', {a}, 'names', {names}, 'p_bar', {p_bar}, 'c', {c}, 'delta', {delta}, ...
                'nan_c', {index_nan_c}, 'nan_b', {index_nan_b}, 'nvi', {nvi}, 'betas', {beta}, ...
                'sectors', {sectors}, 'tkr', {tkr}, 'nvi_benchmark', {nvi_benchmark}, 'uncorr_shocks', {uncorr_shocks}, 'Astar', {Astar_ratio_worst});
            
            results_save_worst(end+1) = results_temp;
    end
    
%     if obj_fun_test
%         results_save_worst(end+1).obj_value_rerun = {};
%         [rerun_v, dl_rerun] = simulation_ratio(xstar,p_bar,c(~index_nan_c),w,g,delta,...
%             N,num_sim*100,use_parallel,index_nan_b,index_nan_c, -1, 0, uncorr_shocks);   
%         results_temp.obj_value_rerun = rerun_v;
%         results_save_worst(end) = results_temp;
%     else
%         results_save_worst(end+1) = results_temp;
%     end
    
   
    %Second run, intial point endpoint of first
    if rerun
        file_name = char(strcat('knitro_rerun.opt'));
        init = xstar;
        % init(init == 0) = 1e-5;
        [xstar,fval,exitflag,output,lambda,grad,t1] = knitro_helper(init,...
            Aineq,bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, ...
            p_bar, c, w, g ,delta, N, num_sim, use_parallel, index_nan_b, ...
            index_nan_c, uncorr_shocks);

        results_save_worst2(end+1) = struct('qt', {q}, 'options', {fileread(file_name)}, ...
            'exit', {exitflag}, 'times', {t1}, 'output', {x000000000star}, ...
            'options_file', {file_name}, 'output_value', {output}, 'obj_value', {fval}, ...
            'w', {w}, 'a', {a}, 'names', {names}, 'p_bar', {p_bar}, 'c', {c}, 'delta', {delta}, ...
            'nan_c', {index_nan_c}, 'nan_b', {index_nan_b}, 'nvi', {nvi}, 'betas', {beta}, ...
            'sectors', {sectors}, 'tkr', {tkr}, 'nvi_benchmark', {nvi_benchmark}, 'uncorr_shocks', {uncorr_shocks});
    end
    
    %Multistart
    if multistart
    file_name = char(strcat('knitro_multistart.opt'));
        sim_test = [1,NaN;10,NaN;100,NaN;1000,NaN;10000,NaN;100000,NaN];
        for j = (3:3)
            num_sim=sim_test(j,1);
            [xstar,fval,exitflag,output,lambda,grad,t1] = knitro_helper(x0,Aineq,...
                bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, p_bar, c, ...
                w, g ,delta, N, num_sim, use_parallel, index_nan_b, index_nan_c, ...
                uncorr_shocks);
            Astar_ratio_worst = reshape(xstar(1:N^2, 1),N,N);
            sim_test(j,2)=fval;
            results_save_worst_multi = struct('qt', {q}, 'options', {fileread(file_name)}, ...
                'exit', {exitflag}, 'times', {t1}, 'output', {xstar}, ...
                'options_file', {file_name}, 'output_value', {output}, 'obj_value', {fval}, ...
                'w', {w}, 'a', {a}, 'names', {names}, 'p_bar', {p_bar}, 'c', {c}, 'delta', {delta}, ...
                'nan_c', {index_nan_c}, 'nan_b', {index_nan_b}, 'nvi', {nvi}, 'betas', {beta}, ...
                'sectors', {sectors}, 'tkr', {tkr}, 'nvi_benchmark', {nvi_benchmark}, 'uncorr_shocks', {uncorr_shocks}, 'Astar', {Astar_ratio_worst});
        end 
    end
    
    %Different numbers of uncorrelated underlying shocks    
    if uncorr_shock_test
    for uncorr_shocks_iter=1:5:N
        file_name = char(strcat('knitro_benchmark.opt'));
        [xstar,fval,exitflag,output,lambda,grad,t1] = knitro_helper(x0,Aineq,...
            bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, p_bar, ...
            c, w, g ,delta, N, num_sim, use_parallel, index_nan_b, index_nan_c, ...
            uncorr_shocks_iter);

        results_save_worst_uncorr_iter(end+1) = struct('qt', {q}, 'options', {fileread(file_name)}, ...
            'exit', {exitflag}, 'times', {t1}, 'output', {xstar}, ...
            'options_file', {file_name}, 'output_value', {output}, 'obj_value', {fval}, ...
            'w', {w}, 'a', {a}, 'names', {names}, 'p_bar', {p_bar}, 'c', {c}, 'delta', {delta}, ...
            'nan_c', {index_nan_c}, 'nan_b', {index_nan_b}, 'nvi', {nvi}, 'betas', {beta}, ...
            'sectors', {sectors}, 'tkr', {tkr}, 'nvi_benchmark', {nvi_benchmark}, 'uncorr_shocks', {uncorr_shocks_iter});
    end


    [xstar,fval,exitflag,output,lambda,grad,t1] = knitro_helper(x0,Aineq,bineq,Aeq,beq,lb,ub,extendedFeatures, options, file_name, p_bar, c, w, g ,delta, N, num_sim, use_parallel, index_nan_b, index_nan_c, N);
    results_save_worst_uncorr_iter(end+1) = struct('qt', {q}, 'options', {fileread(file_name)}, ...
        'exit', {exitflag}, 'times', {t1}, 'output', {xstar}, ...
        'options_file', {file_name}, 'output_value', {output}, 'obj_value', {fval}, ...
        'w', {w}, 'a', {a}, 'names', {names}, 'p_bar', {p_bar}, 'c', {c}, 'delta', {delta}, ...
        'nan_c', {index_nan_c}, 'nan_b', {index_nan_b}, 'nvi', {nvi}, 'betas', {beta}, ...
        'sectors', {sectors}, 'tkr', {tkr}, 'nvi_benchmark', {nvi_benchmark}, 'uncorr_shocks', {N});
    end


    %Different tests of numberic gradient...
    % gradcheck_num = 10000;
    % if gradcheck
    %     disp("Starting Gradient Checks...")
    %     tester = 0;
    %     stopper = 0;
    %     disp("Starting standard. Forward differences, optimal (worst) network...")
    %     j = 1;
    %     while stopper == 0
    %         try
    %         tester = tester + 1;
    %         tic
    %         obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,gradcheck_num,use_parallel,index_nan_b,index_nan_c, -1); 
    %         [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),results_save_worst(1).output,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_gradcheck_forward.opt');%,'knitro.opt'
    %         t1 = toc;
    %         stopper = 1;
    %         catch err
    %             disp(tester)
    %         end
    %     end
    % 
    %     tester = 0;
    %     stopper = 0;
    %     disp("Starting standard. Forward differences, optimal (worst) network, with more simulations (50000)...")
    %     j = 1;
    %     while stopper == 0
    %         try
    %         tester = tester + 1;
    %         tic
    %         obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,50000,use_parallel,index_nan_b,index_nan_c, -1); 
    %         [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),results_save_worst(1).output,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_gradcheck_forward.opt');%,'knitro.opt'
    %         t1 = toc;
    %         stopper = 1;
    %         catch err
    %             disp(tester)
    %         end
    %     end
    %     
    %     
    %     tester = 0;
    %     stopper = 0;
    %     j = 1;
    %     disp("Reaching feasible but non-optimal network...")
    %     while stopper == 0
    %         try
    %         tester = tester + 1;
    %         tic
    %         obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,gradcheck_num,use_parallel,index_nan_b,index_nan_c, -1); 
    %         [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),x0,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_gradcheck_stop20.opt');%,'knitro.opt'
    %         t1 = toc;
    %         stopper = 1;
    %         catch err
    %             disp(tester)
    %         end
    %     end
    %     
    %     
    %     tester = 0;
    %     stopper = 0;
    %     j = 1;
    %     disp("Gradient check over non-optimal network (forward)...")
    %     while stopper == 0
    %         try
    %         tester = tester + 1;
    %         tic
    %         obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,gradcheck_num,use_parallel,index_nan_b,index_nan_c, -1); 
    %         [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),xstar,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_gradcheck_forward.opt');%,'knitro.opt'
    %         t1 = toc;
    %         stopper = 1;
    %         catch err
    %             disp(tester)
    %         end
    %     end
    %     
    %     disp("Starting optimal network with central differences...")
    %     tester = 0;
    %     stopper = 0;
    %     j = 1;
    %     while stopper == 0
    %         try
    %         tester = tester + 1;
    %         tic
    %         obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,gradcheck_num,use_parallel,index_nan_b,index_nan_c, -1); 
    %         
    %         [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),results_save_worst(1).output,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_gradcheck_central.opt');%,'knitro.opt'
    %         
    %         t1 = toc;
    %         stopper = 1;
    %         catch err
    %             disp(tester)
    %         end
    %     end
    % end

    
    
    save([to_save_file '_intermed'])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_network == 1
results_to_use = results_save_worst;
%results_to_use = results_save_worst_multi;
num_qt = size(results_to_use, 2);
tile = ceil(sqrt(num_qt));
for to_plot=1:num_qt
%     if results_to_use(to_plot).uncorr_shocks == 1
    nodes_in_graph = N;
    cutoff = 1e-2; % do not plot connections smaller than cutoff
    % get relative liabilities matrix
    Astar_ratio_worst = reshape(results_to_use(to_plot).output(1:N^2, 1),N,N);
    % Astar_ratio_best = reshape(results_save_best(to_plot).output(1:N^2, 1),N,N);

    % construct normalizing constant to convert to same units across all plots
    Astar_ratio_worst_abs = Astar_ratio_worst .* repmat(results_to_use(to_plot).p_bar', 1, N);
    % Astar_ratio_best_abs = Astar_ratio_best .* repmat(results_save_worst(to_plot).p_bar', 1, N);
    max_link = max(max([Astar_ratio_worst_abs]));

    % worst network plot
    f1=figure(to_plot);
    f1.Color=[1 1 1];
    f1.Position = [0 0 1000 1000];

    Astar_ratio_worst_abs(Astar_ratio_worst_abs<cutoff)=0; % ignore small links
    Astar_ratio_worst_abs = Astar_ratio_worst_abs(1:nodes_in_graph, 1:nodes_in_graph);
    [temp, sector_sort] = sort(results_to_use(to_plot).sectors(1:nodes_in_graph));
    Astar_ratio_worst_abs = Astar_ratio_worst_abs(sector_sort, sector_sort);
    sectors_unique = unique(results_to_use(to_plot).sectors(sector_sort));
    
    Gstar = digraph(Astar_ratio_worst_abs); % create directed graph
    p1=plot(Gstar,'Layout','circle');% options for layout include: 'layered','circle','force','subspace','force3','subspace3'
    Gstar.Edges.LWidths =7*Gstar.Edges.Weight/max_link; % make edge width proportional to link size in dollars (normalized by max_link)
    p1.LineWidth = Gstar.Edges.LWidths;
    p1.NodeLabel = results_to_use(to_plot).tkr(sector_sort);
    %p1.MarkerSize = 30*results_save_worst(to_plot).a/max(results_save_worst(to_plot).a); % make node size proportional to assets in dollars (normalized by max_link)
    p1.MarkerSize = 30*results_to_use(to_plot).w(sector_sort)/max(results_to_use(to_plot).w); % make node size proportional to equity capital in dollars (normalized by max_link)
    p1.NodeColor = 'r';

    % delete arrows and node labels, which are created manually below
    p1.ArrowSize=0;
    nl = p1.NodeLabel;
    p1.NodeLabel = '';

    % get position of nodes
    xd = get(p1, 'XData');
    yd = get(p1, 'YData');

    % create x-y coordinates for node labels
    xd_text = xd;
    yd_text = yd;
    xd_text(xd_text>0)=xd_text(xd_text>0)*1.05;
    xd_text(xd_text<0)=xd_text(xd_text<0)*1.20;
    yd_text(yd_text>0)=yd_text(yd_text>0)*1.15;
    yd_text(yd_text<0)=yd_text(yd_text<0)*1.15;

    sectors_to_use = results_to_use(to_plot).sectors(sector_sort);
    for j=1:size(sectors_unique,2)
       nl(end+1) = sectors_unique(j);
       ix = strcmp(sectors_to_use, sectors_unique(j));
       x_coord = xd_text(ix);
       x_coord = x_coord(1)*1.2;
       y_coord = yd_text(ix);
       y_coord = y_coord(1)*1.2;
       xd_text(end+1) = x_coord;
       yd_text(end+1) = y_coord;
    end

    % plot node labels
    text(xd_text,yd_text, nl, 'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')

    % get x-y coordinates of edges (first column is origin, second column is
    % end)
    [E1,E2]=find(Astar_ratio_worst_abs'>0);
    arrow_XData =[xd(E1)' xd(E2)'];
    arrow_YData =[yd(E1)' yd(E2)'];

    % reduce white space around figure
    ax = f1.CurrentAxes;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height]*0.95;

    %get axes and drawing area limits
    ax_lim = ax.Position;
    ax_xlim = [ax_lim(1) ax_lim(3)+ax_lim(1)];
    ax_ylim = [ax_lim(2) ax_lim(4)+ax_lim(2)];
    fig_xlim = xlim(ax);
    fig_ylim = ylim(ax);

    % convert coordinates for edges to units that matlab understands (between 0 and 1)
    Xnorm = (arrow_XData-fig_xlim(1))*diff(ax_xlim)/diff(fig_xlim)+ax_xlim(1);
    Ynorm = (arrow_YData-fig_ylim(1))*diff(ax_ylim)/diff(fig_ylim)+ax_ylim(1);

    % plot arrows
    lambda = 0.65*ones(size(Xnorm,1),1); % vector of arrow locations within each edge; 0.5 plots the arrow in the middle of the edge
    for i=1:size(Xnorm,1)
        annotation('arrow','Position',[Xnorm(i,2) Ynorm(i,2) (lambda(i)-1)*diff(Xnorm(i,:)) (lambda(i)-1)*diff(Ynorm(i,:))],'LineStyle','none','Color','k','HeadLength',max(p1.LineWidth(i),3),'HeadWidth',max(p1.LineWidth(i),3))
    end
    axis off


    % save with minimal margin
    f1.PaperPositionMode = 'auto';
    fig_pos = f1.PaperPosition;
    f1.PaperSize = [fig_pos(3) fig_pos(4)];

    print(f1, ['Figures/network_schem_worst_' spec_id],'-dpdf')
    end
% end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Plot tests of uncorrelated shock numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beta
% 192 .7576203
% 193 .7318331
% 194 .7506473
% 195 .7025943
% 196 .7160696
% 197 .7006595
% 198 .6984555
% 199 .6883392
% 200 .70696151

if plot_uncorrs
    res_uncorr_test = nan(q_end-q_start+1, 50);

    for iter=1:size(results_save_worst_uncorr_iter, 2)
        q = results_save_worst_uncorr_iter(iter).qt
        res_uncorr_test(q-q_start+1, results_save_worst_uncorr_iter(iter).uncorr_shocks) = ...
            100*((-1 * results_save_worst_uncorr_iter(iter).obj_value)-1);        
    end

    labels = 1:50
    labels = labels(~any(isnan(res_uncorr_test)));
    res_uncorr_test(:, any(isnan(res_uncorr_test))) = [];
    qt = [q_start:q_end]';
    f = figure(10)
    plot(qt, res_uncorr_test)
    legend(cellstr(string(labels)))
     xticks(q_start:q_end)
    ylabel('%')
    set(get(gca, 'ylabel'), 'rotation', 0)
    set(gcf, 'Color', [1, 1, 1])
    print(f,['Figures/nvi_compare_uncorr_shocks_all_' spec_id],'-dpdf')
    
    f = figure(11)
    plot(qt, res_uncorr_test(:, [1, 3, 6, 9, 11]))
    legend(cellstr(string(labels([1, 3, 6, 9, 11]))))
     xticks(q_start:q_end)
    ylabel('%')
    set(get(gca, 'ylabel'), 'rotation', 0)
    set(gcf, 'Color', [1, 1, 1])
    print(f,['Figures/nvi_compare_uncorr_shocks_' spec_id],'-dpdf')
        
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Plot obj fun test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_obj_fun_test
    sim_results = nan(q_end-q_start+1, 1);
    sim_results_obj_rerun = nan(q_end-q_start+1, 1);
    for q=q_start:q_end
        sim_results_obj_rerun(q-q_start+1) = 100*((-1 * results_save_worst(q-q_start+1).obj_value_rerun) - 1);
        sim_results(q-q_start+1) = 100*((-1 * results_save_worst(q-q_start+1).obj_value) - 1);
    end
    
    qt = [q_start:q_end]';
    f = figure(9)
    plot(qt, sim_results, qt, sim_results_obj_rerun)
    legend({'Simulation Obj Value', 'Simulation Obj Value, Re-run (Num sim x 100)'}, 'Location', 'southoutside')
    xticks(q_start:q_end)
    ylabel('%')
    set(get(gca, 'ylabel'), 'rotation', 0)
    set(gcf, 'Color', [1, 1, 1])
    print(f,['Figures/obj_fun_test_' spec_id],'-dpdf')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Plot obj fun and NVIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beta_max -- gets hard-coded in, for now. Will be applied to post-sim
% network (otherwise values explode)
% 192 .7576203
% 193 .7318331
% 194 .7506473
% 195 .7025943
% 196 .7160696
% 197 .7006595
% 198 .6984555
% 199 .6883392
% 200 .70696151

if plot_nvi
    nvi = nan(q_end-q_start+1, 1);
    nvi_after_benchbeta = nan(q_end-q_start+1, 1);
    nvi_after_realbeta = nan(q_end-q_start+1, 1);
    sim_results = nan(q_end-q_start+1, 1);
    sim_results2 = nan(q_end-q_start+1, 1);
    sim_results_obj_rerun = nan(q_end-q_start+1, 1);
    sim_results_multi = nan(q_end-q_start+1, 1);
    beta_pluses = nan(q_end-q_start+1, 1);
    beta_pluses_bench = nan(q_end-q_start+1, 1);
    for q=q_start:q_end
            if q == 192
               beta_pluses_bench(q-q_start+1) = .7576203;
            elseif q == 193
               beta_pluses_bench(q-q_start+1) = .7318331;
            elseif q == 194
                beta_pluses_bench(q-q_start+1) = .7506473;
            elseif q == 195
                beta_pluses_bench(q-q_start+1) = .7025943;
            elseif q == 196
                beta_pluses_bench(q-q_start+1) = .7160696;
            elseif q == 197
                beta_pluses_bench(q-q_start+1) = .7006595;
            elseif q == 198
                beta_pluses_bench(q-q_start+1) = .6984555;
            elseif q == 199 
                beta_pluses_bench(q-q_start+1) = .6883392;
             elseif q == 200 
                beta_pluses_bench(q-q_start+1) = .70696151;
            end
            nvi(q-q_start+1) = max(results_save_worst(q-q_start+1).nvi_benchmark);
            sim_results(q-q_start+1) = 100*((-1 * results_save_worst(q-q_start+1).obj_value) - 1);
%             sim_results2(q-q_start+1) = 100*((-1 * results_save_worst2(q-q_start+1).obj_value) - 1);
            if multistart
                sim_results_multi(q-q_start+1) = 100*((-1 * results_save_worst_multi(q-q_start+1).obj_value) - 1);
                Astar_ratio_worst_temp = reshape(results_save_worst_multi(q-q_start+1).output(1:N^2, 1),N,N);
            
            else
                Astar_ratio_worst_temp = reshape(results_save_worst(q-q_start+1).output(1:N^2, 1),N,N);
            end
            betas = sum(Astar_ratio_worst_temp, 2);
            beta_plus = max(betas);
            beta_pluses(q-q_start+1) = beta_plus;
            
            c_actual = results_save_worst(q-q_start+1).c;
            c_sim = results_save_worst(q-q_start+1).output(N^2+sum(results_save_worst(q-q_start+1).nan_b)+1:end);
            c_actual(results_save_worst(q-q_start+1).nan_c) = c_sim;
            delta_a = results_save_worst(q-q_start+1).a .* results_save_worst(q-q_start+1).delta;
            nvi_after_benchbeta(q-q_start+1) = 100 * sum(delta_a)/((1-(1+g)*beta_pluses_bench(q-q_start+1))*sum(results_save_worst(q-q_start+1).a));
            nvi_after_realbeta(q-q_start+1) = 100 * sum(delta_a)/((1-(1+g)*beta_pluses(q-q_start+1))*sum(results_save_worst(q-q_start+1).a));

    end
    qt = [q_start:q_end]';
    f = figure(10)
    plot(qt, nvi, qt, sim_results, qt, sim_results_multi, qt, nvi_after_benchbeta, qt, nvi_after_realbeta)
    legend({'Benchmark NVI', 'Simulation', 'Simulation with Multistart', 'NVI with Post-Sim Network Benchmark Beta', 'NVI with Post-Sim Network Real Beta'}, 'Location', 'southoutside')
    xticks(q_start:q_end)
    ylabel('%')
    set(get(gca, 'ylabel'), 'rotation', 0)
    set(gcf, 'Color', [1, 1, 1])
    print(f,['Figures/nvi_compare_' spec_id],'-dpdf')

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output final results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -EXCEPT results_save_worst results_save_worst2 results_save_worst ...
    results_save_worst_multi results_save_best q_start q_end ...
    num_sim use_parallel N g to_save_file uncorr_shocks spec_id ...
    multistart num_sim run_analytics_tests run_nvi_tightness


save(to_save_file)
diary off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Some Analytics testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_analytics_tests
    %Assuming first node is the connecting one, here
    w = results_save_worst_multi.w;
    c = results_save_worst_multi.c;
    c(results_save_worst_multi.nan_c) = results_save_worst_multi.output(N^2 + sum(results_save_worst_multi.nan_b)+1:end);

    a12 = results_save_worst_multi.Astar(1, 2);

    d1 = 1-w(1)/c(1);
    d2 = 1 - w(2)/c(2);

    prob1_not_2 = w(2)^2;
    prob1_not_2 = prob1_not_2 / (2*a12*c(1)*c(2));

    prob_both_cause1 = 2*(a12*c(1)*w(2) - w(2)*w(2) - w(1)*w(2)*a12) + w(2)^2;
    prob_both_cause1 = prob_both_cause1 / (2*a12*c(1)*c(2));

    prob_both_indep = c(2)*(c(1)-w(1)) - w(2)*(c(1)-w(1));
    prob_both_indep = prob_both_indep / (c(1)*c(2));

    prob_both_1small = w(2)*(2*c(2)-w(2));
    prob_both_1small = prob_both_1small/(2*a12*c(1)*c(2));

    prob_both_1big = (c(1) - (w(2)/a12 + w(1)))/c(1);
    prob_both = prob_both_1small + prob_both_1big;

    test_size = 100000000;
    sim_test = random('Uniform',0,1,[test_size,N]).*repmat(c,test_size,1);

    test_prob1_not_2 = sum(sim_test(:, 1) > w(1) & a12*(sim_test(:, 1) - w(1)) + sim_test(:, 2) < w(2)) / test_size;

    %% 
    %Exp of x1 - w1, given that node 1 defaults but node 2 does not

    exp_1 = w(2)^3;
    exp_1 = exp_1 / (6 * a12^2 * c(1) * c(2));
    % exp_1 = exp_1 / prob1_not_2;

    points = sim_test(:, 1) > w(1) & a12*(sim_test(:, 1) - w(1)) + sim_test(:, 2) < w(2);
    exp_test = mean(sim_test(points, 1)) - w(1);

    %%
    %Exp of x1 - w1, given that both nodes default -- and node 1's shock
    %doesn't automatically with node 2 out

    exp_2 = w(2)^2 * (3*c(2) - w(2));
    exp_2 = exp_2 / (6 * a12 ^2 * c(1) * c(2));
    % exp_2 = exp_2 / (prob_both_1small);

    points = sim_test(:, 1) > w(1) & a12*(sim_test(:, 1) - w(1)) +  sim_test(:, 2) > w(2) & sim_test(:, 1) < (w(2)/a12+w(1));
    exp_test = mean(sim_test(points, 1)) - w(1);


    %%
    %Exp of x1-w1, given that node 1's shock is big enough to automatically wip
    %out node 2

    exp_3 = a12^2*(c(1)-w(1))^2 - w(2)^2;
    exp_3 = exp_3 / (2*a12^2 * c(1));
    % exp_3 = exp_3 / (prob_both_1big);

    points = sim_test(:, 1) > (w(2)/a12+w(1));
    exp_test = mean(sim_test(points, 1)) - w(1);


    Lbar_true = c(2)/2 + w(1)*d1 + (1-d1)*w(1)/2 + exp_1 + (1+a12)*(exp_2+exp_3);

    L0 = c(1)/2 + c(2)/2;

    ratio_true = Lbar_true/L0;

    nvi_here = (c(1) * d1 + c(2)*d2)/(c(1) + c(2));
    nvi_here = nvi_here / (1-a12);
    nvi_here = nvi_here + 1;
end

if run_nvi_tightness

    [rerun_v, dl_rerun] = simulation_ratio(results_save_worst_multi.output,results_save_worst_multi.p_bar, ...
        results_save_worst_multi.c(~results_save_worst_multi.nan_c),results_save_worst_multi.w,g,results_save_worst_multi.delta,...
        N,100*num_sim,use_parallel,results_save_worst_multi.nan_b,results_save_worst_multi.nan_c, -1, 0, uncorr_shocks, 'uniform');   

    
end





