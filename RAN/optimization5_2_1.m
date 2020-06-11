
clear
use_parallel = 0; 
pool_size = 90;
if use_parallel
    p = gcp('nocreate'); 
    if isempty(p)
         % if pool not open, create one
         mypool=parpool(pool_size,'AttachedFiles',{'ui.m','systemeq.m'});
    end
end

%% Set up optimization

% Network primitives
p_bar = [55 55 140 55 55]; % total liabilities
c = [50 50 150 50 50]; % outside assets
b = [55 55 100 55 55]; % outside liabilities
w = [5 5 10 5 5]; % net worth
g = 0.05; % bankruptcy cost

% Other network variables
a = w+p_bar; % total assets
d = a-c; % inside assets
f = p_bar-b; % inside liabilities
N = length(c); % number of nodes

% Number of trials for simulation
num_sim = 5;

% Initial guess 
A0 = [ 0 0 0 0 0; 0 0 0 0 0; 1/14 1/14 0 1/14 1/14; 0 0 0 0 0; 0 0 0 0 0];
A0 = A0(:);

% linear inequality constraints 
Aineq = repmat(eye(N),1,N); bineq = ones(N,1); % sum of rows of A cannot exceed 1

% linear equality constraints 
Aeq1=kron(eye(N),p_bar); beq1=d';  % for each i, sum( p_bar * ith column of A) = inside assets
Aeq2=repmat(diag(p_bar),1,N); beq2 =f';  %for each i, p_bar_i*sum( ith row of A) = inside liabilities
Aeq3=full(sparse(1:5,sub2ind([N,N],1:N,1:N),1,N,N^2));beq3=zeros(N,1); % diagonal elements of A must be zero

Aeq = [Aeq1;Aeq2;Aeq3];
beq = [beq1;beq2;beq3]; 

% upper and lower bounds: elements of A between 0 and 1
lb = zeros(size(A0));
ub = ones(size(A0));

% non-linear constraints
nonlcon = @(A) sums(A,N);

%% check if user-provided gradient matches with numerical one computed by Matlab

fmincon(@(A) simulation(A,p_bar,c,w,g,N,num_sim,use_parallel),A0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,optimoptions('fmincon','Display','iter','MaxIter',1,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',true,'FinDiffType','central'));

%%   Use knitro

%options = optimset('Display','iter','TolX',1e-2,'TolFun',1e-2,'DiffMinChange',0.01);
% [Astar,fval,exitflag,output,lambda]  = knitromatlab( @(A) simulation(A,p_bar,c,w,g,N,num_sim,use_parallel),A0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,[],options,'knitro.opt');
% Astar = reshape(Astar,N,N);
% 

%% Use matlab global search 
% options = optimoptions(@fmincon,'Display','iter-detailed','StepTolerance',1e-3,'FunctionTolerance',0.01,'ConstraintTolerance',0.001,'OptimalityTolerance',0.01,...
%     'FiniteDifferenceStepSize',0.001,'MaxFunctionEvaluations', 10000, 'MaxIterations', 500,...
%     'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,'CheckGradients',false,...%'FinDiffType','central',...
%     'Algorithm','sqp','ScaleProblem','obj-and-constr');%,'TypicalX',x0+10^(-7));
% problem = createOptimProblem('fmincon','objective',...
%     @(A) simulation(A,p_bar,c,w,g,N,num_sim,use_parallel),'x0',A0,'Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nonlcon(x),'options',options);
% gs = GlobalSearch('Display','iter','FunctionTolerance',0.1,'XTolerance',0.1);
% 
% [Astar,fval,exitflag,output,solutions] = run(gs,problem);
% Astar = reshape(Astar,N,N);

%% Use matlab multistart 

options = optimoptions(@fmincon,'Display','iter-detailed','StepTolerance',1e-3,'FunctionTolerance',0.01,'ConstraintTolerance',0.001,'OptimalityTolerance',0.01,...
    'FiniteDifferenceStepSize',0.001,'MaxFunctionEvaluations', 10000, 'MaxIterations', 500,...
    'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,'CheckGradients',false,...%'FinDiffType','central',...
    'Algorithm','sqp','ScaleProblem','obj-and-constr');%,'TypicalX',x0+10^(-7));
problem = createOptimProblem('fmincon','objective',...
    @(A) simulation(A,p_bar,c,w,g,N,num_sim,use_parallel),'x0',A0,'Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nonlcon(x),'options',options);
ms = MultiStart('Display','iter','FunctionTolerance',0.1,'XTolerance',0.1,'StartPointsToRun','bounds-ineqs');

[Astar,fval,exitflag,output,solutions]  = run(ms,problem,5);
Astar = reshape(Astar,N,N);

%% Use matlab pattern search 

% options = optimset('Display','iter','TolX',1e-4,'TolFun',1e-2,'DiffMinChange',0.001);
% [Astar,fval,exitflag,output] = patternsearch(@(A) simulation(A,p_bar,c,w,g,N,num_sim,use_parallel),A0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
% Astar = reshape(Astar,N,N);

  
% [x2,FVAL2,EXITFLAG2,OUTPUT2,LAMBDA2]= fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optimset('Display','iter'));
% x2 = reshape(x2,sqrt(length(x2)),sqrt(length(x2)))

%% Close parallel pool, save

% delete(gcp('nocreate'))
save('sim_results')