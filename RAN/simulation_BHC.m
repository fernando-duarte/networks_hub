%TEST TEST THIS IS A TEST
clear
addpath(genpath('Toolboxes'))
diary on
% plotlysetup('fernando.duarte', 'ReAHRvXaOy06uPcsIZKJ')
% saveplotlycredentials('fernando.duarte', 'ReAHRvXaOy06uPcsIZKJ')

%% Set up optimization

sheets ={'Broker Dealers', 'Insurance','Other','REITs','BHCs'};%,'Approx Aggregates'
nsheets = length(sheets);
load('network_data_196')
% if isunix                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
%     load('network_data')
% else
%     num=cell(nsheets,1);txt=cell(nsheets,1);raw=cell(nsheets,1);headers=cell(nsheets,1);names=cell(nsheets,1);
%     p_bar=cell(nsheets,1);c=cell(nsheets,1);b=cell(nsheets,1);w=cell(nsheets,1);delta=cell(nsheets,1);
%     a=cell(nsheets,1);d=cell(nsheets,1);f=cell(nsheets,1);N_tot=cell(nsheets,1);type=cell(nsheets,1);
%     qt=cell(nsheets, 1);
%     
%     for i=1:nsheets
%         [num{i},txt{i},raw{i}] = xlsread('node_stats_for_simulation.xlsx',sheets{i});
%         headers{i} = raw{i}(1,:);
%         names{i} = raw{i}(2:end,1);
%         
%         % Network primitives
% %         p_bar{i} = vertcat(raw{i}{2:end,strcmpi('p_bar',headers{i})}); % total liabilities
%         c{i} = vertcat(raw{i}{2:end,strcmpi('c',headers{i})}) ; % outside assets
%         b{i} = vertcat(raw{i}{2:end,strcmpi('b',headers{i})}); % outside liabilities
%         w{i} =vertcat(raw{i}{2:end,strcmpi('w',headers{i})}); % net worth
%         delta{i} =vertcat(raw{i}{2:end,strcmpi('delta',headers{i})}); % probability of default
%         a{i} = vertcat(raw{i}{2:end,strcmpi('assets',headers{i})});  % total assets
%         qt{i} = vertcat(raw{i}{2:end,strcmpi('qt_dt',headers{i})});  % total assets
%         
%         % Other network variables
% %         a{i} = w{i}+p_bar{i}; % total assets
%         p_bar{i} =a{i} - w{i}; % total liabilities
%         d{i}=  a{i}-c{i};% inside assets
%         f{i} = p_bar{i}-b{i};% inside liabilities
%         N_tot{i} = length(c{i}); % number of nodes
%         type{i} = repmat(i,N_tot{i},1); % type of firm
%         save('network_data')
%     end
% end

scale = 1e6;
p_bar = cell2mat(p_bar)/scale;
c = cell2mat(c)/scale;
b = cell2mat(b)/scale;
w = cell2mat(w)/scale;
delta = cell2mat(delta);
a = cell2mat(a)/scale;
d = cell2mat(d)/scale;
f = cell2mat(f)/scale;
type = cell2mat(type);
names = vertcat(names{:});
 
% drop if missing total assets
keep = ~isnan(a);
a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep);
names = names(keep);

% drop if net worth is negative
keep = w>0;
a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep);
names = names(keep);

% drop if any variable is missing
% keep = ~any(isnan([a c b w delta d f p_bar]),2);
% a = a(keep);c = c(keep);b= b(keep);w = w(keep);delta= delta(keep);d = d(keep);f = f(keep);type = type(keep);p_bar =p_bar(keep);
% names = names(keep);

% sort by total assets
[a,ix] = sortrows(a,-1);
c = c(ix);b= b(ix);w = w(ix);delta= delta(ix);d = d(ix);f = f(ix);type = type(ix);p_bar =p_bar(ix);
names = names(ix);

%% open parallel pool
use_parallel = 1; 
pool_size = 30;
if use_parallel
    p = gcp('nocreate'); 
    if isempty(p)
         % if pool not open, create one
         mypool=parpool(pool_size,'AttachedFiles',{'ui.m','systemeq.m','sums.m'});
    end
end
%%
% Number of trials for simulation
num_sim = 1000;

% number of nodes to keep (top N by assets)
N = 10;%length(a);%795;

% bankruptcy cost
g = 0.05; 

a = a(1:N)';c = c(1:N)';b= b(1:N)';w = w(1:N)';delta= delta(1:N)';d = d(1:N)';f = f(1:N)';type = type(1:N)';p_bar =p_bar(1:N)';
names = names(1:N)';

% find missing values for b and c
index_nan_b = isnan(b);
index_nan_c = isnan(c);
index_nan  = index_nan_b | index_nan_c;

%% glasserman 
% 
% y=10;
% p_bar = [55+y 55+y 140 55+y 55+y]; % tot liab
% A0 = [0 y/p_bar(1) 0 0 0; 0 0 0 0 y/p_bar(1); 1/14 1/14 0 1/14 1/14; y/p_bar(4) 0 0 0 0; 0 0 0 y/p_bar(5) 0];
% A0 = A0(:);
% 
% c = [50 50 150 50 50]; % outside assets
% w = [5 5 10 5 5]; % net worth
% g = 0.05;
% 
% % Network primitives
% b = [55 55 100 55 55]; % outside liabilities
% a = w+p_bar; % total assets
% d=  a-c;% inside assets
% f = p_bar-b;% inside liabilities
% N = length(c); % number of nodes
% 
%%

% linear inequality constraints 
Aineq = [];bineq=[];%Aineq=repmat(eye(N),1,N); bineq = (f./p_bar)'; % sum of columns of A cannot exceed ratio of inside to total liabilities

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

% Aeq1=[];beq1=[];
%check
% [Aeq1*A0 beq1]-[sum(reshape(A0,N,N).*repmat(p_bar',1,N),1); a-c]'
% [Aeq2*A0 beq2]-[sum(reshape(A0,N,N),2).*p_bar' p_bar'-b']
% [w-c+p_bar; sum(reshape(A0,N,N).*repmat(p_bar',1,N),1)]


% upper and lower bounds:
% 1) elements of A between 0 and ratio of inside to total liabilities
% 2) elements of b0 between zero+1e-5 and p_bar - 1e-5
% 3) elements of c0 between w+1e-5 and a - 1e-5
lb_A=sparse(N^2,1);   ub_A = repmat((f./p_bar)',N,1); ub_A(isnan(ub_A))=1;
lb_b=sparse(sum(index_nan_b),1)+1e-5;   ub_b = ones(sum(index_nan_b),1).*p_bar(index_nan_b)'- 1e-5;
lb_c=ones(sum(index_nan_c),1).*w(index_nan_c)'+1e-5;   ub_c = ones(sum(index_nan_c),1).*a(index_nan_c)'-1e-5;

lb = [lb_A;lb_b;lb_c];
ub = [ub_A;ub_b;ub_c];

%check
%[A0 ub]

% non-linear constraints (complementary constraints implemented as
% non-linear function)
nonlcon = @(x) sums(x,N,sum(index_nan_b)+sum(index_nan_c));

% complementarity constraints
upper_triang = ~tril(ones(N),0);
lower_triang = ~triu(ones(N),0);
[ix_row,ix_column]=find(upper_triang);
ix_transpose = fliplr([ix_row ix_column]);
ix =sortrows([sub2ind([N N],ix_transpose(:,1),ix_transpose(:,2)) find(upper_triang)],2);

extendedFeatures.ccIndexList1 = ix(:,1)';
extendedFeatures.ccIndexList2 = ix(:,2)';

%% initial guess for A

% reshape(Aeq\beq,N,N)
% reshape(pinv(Aeq)*beq,N,N)
% [x0,resnorm,residual,exitflag,output,lambda] = lsqnonneg(Aeq,beq,optimset('TolX',0.001));
% lb_violation = x0<=lb;
% ub_violation = x0>=ub;
% 
% x0(lb_violation)=lb(lb_violation)+(ub(lb_violation)-lb(lb_violation))/2;
% x0(ub_violation)=lb(ub_violation)+(ub(ub_violation)-lb(ub_violation))/2;
% A0 = reshape(x0(1:N^2,N,N));  
% b0 = x0(N^2+1:N^2+sum(index_nan_b));  
% c0 = x0(N^2+sum(index_nan_b)+1:end);  

% upper_triang = ~tril(ones(size(A0)),1);
% length_x0 = N^2+sum(index_nan_b)+sum(index_nan_c);
% x0 = sparse([find(upper_triang);(N^2+1:N^2+sum(index_nan_b))';(N^2+sum(index_nan_b)+1:length(x0))'],1,[A0(upper_triang);b0;c0],length(length_x0),1,floor(length(length_x0)/2));

x0 = lb+(ub-lb)/2;


% plot initial guess
% G = digraph(reshape(A0,N,N)); % digraph,names(1:N)
% p1=plot(G,'EdgeLabel',round(G.Edges.Weight,2));
% fig = fig2plotly(gcf,'offline',true,'filename','offline-graph');

% G.Edges.Weight = round(G.Edges.Weight,2);
% G.Edges.LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);
% h.LineWidth = G.Edges.LWidths;

%% plot
% Nplot=min(150,N);
% G0 = digraph(A0(1:Nplot,1:Nplot),names(1:Nplot)); % digraph
% p2=plot(G0,'Layout','force');%'layered','circle','force','subspace','force3','subspace3'
% G0.Edges.LWidths =2*G0.Edges.Weight/max(G0.Edges.Weight);
% p2.LineWidth = G0.Edges.LWidths;
% p2.MarkerSize = 10*w(1:Nplot)/max(w(1:Nplot));
% 
% p2.NodeColor = repmat([1 0 0],Nplot,1)
% % p2.NodeColor(type'==5,:) = repmat([1 0 0],sum(type==5),1);
% p2.NodeColor(type(1:Nplot)'==2,:) = repmat([1 1 1],sum(type(1:Nplot)==2),1);
% p2.NodeColor(type(1:Nplot)'==4,:) = repmat([0 1 1],sum(type(1:Nplot)==4),1);
% 
% p2.EdgeAlpha=0.3
% 
% layout(p2,'force','Iterations',30) 
% layout(p2,'subspace','dimension',30) 
% 
% layout(p2,'layered','Direction','right','Sources',[1 4])
%% check if user-provided gradient matches with numerical one computed by Matlab
% options=optimoptions('fmincon','Display','iter','MaxIter',1,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false,'CheckGradients',true,'FinDiffType','central');
% fmincon(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,1,0,index_nan_b,index_nan_c),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
% fmincon(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,1,0,index_nan_b,index_nan_c),x0,[],[],[],[],[],[],[],options);

% 
% options = optimset('FinDiffType','central', 'GradObj','on','GradConstr','off','DerivativeCheck','on','MaxIter',1);
% knitromatlab( @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,1,0,index_nan_b,index_nan_c),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,[],options,'knitro.opt');
% knitromatlab( @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,1,0,index_nan_b,index_nan_c),x0,[],[],[],[],lb,ub,[],extendedFeatures,options);


%%   Use knitro
%,'DiffMinChange',0.005
diary(['diary_' datestr(now) '_' char(string(num_sim)) 'trials_' char(string(N)) 'nodes'])
iterations_to_run = 1;
results_specs = struct('options', {''}, 'exit', zeros(iterations_to_run, 1), 'times', zeros(iterations_to_run, 1), 'output', {''}, 'options_files', {''}, 'output_values', {''}, 'obj_value', zeros(iterations_to_run, 1));
results_specs_ratio = struct('options', {''}, 'exit', zeros(iterations_to_run, 1), 'times', zeros(iterations_to_run, 1), 'output', {''}, 'options_files', {''}, 'output_values', {''}, 'obj_value', zeros(iterations_to_run, 1));

options = optimset('Display','iter-detailed','TolX',1e-4,'TolFun',1e-2, 'TolCon', 1e-3,'GradConstr','off', 'GradObj','on');%,
% options = optimset('Display','iter-detailed', 'GradObj','on','GradConstr','on');


j = 1;
for i=1:1

file_name = char(strcat('knitro_me_', string(i), '.opt'));

% tester = 0;
% while results_specs.times == 0
%     try
%     tester = tester + 1;
%     tic
%     obj = @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c, -1); 
%     [xstar,fval,exitflag,output,lambda,grad]  = knitromatlab( @(x) obj(x),x0,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, file_name);%,'knitro.opt'
%     t1 = toc;
%     results_specs.output{j} = xstar;
%     results_specs.options{j} = fileread(file_name);
%     results_specs.options_files{j} = file_name;
%     results_specs.output_values{j} = output;
%     results_specs.times(j) = t1;    
%     results_specs.obj_value(j) = fval;
%     results_specs.exit(j) = exitflag;
%     catch err
%         disp(tester)
%     end
% end

% tester = 0;
% while results_specs_ratio.times == 0
%     try 
%         tester = tester + 1;
%     tic
%     obj = @(x) simulation_ratio(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c, -1); 
%     [xstar_ratio,fval_ratio,exitflag_ratio,output_ratio,lambda,grad]  = knitromatlab( @(x) obj(x),x0,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, file_name);%,'knitro.opt'
%     t1 = toc;
%     results_specs_ratio.output{j} = xstar_ratio;
%     results_specs_ratio.options{j} = fileread(file_name);
%     results_specs_ratio.options_files{j} = file_name;
%     results_specs_ratio.output_values{j} = output_ratio;
%     results_specs_ratio.times(j) = t1;    
%     results_specs_ratio.obj_value(j) = fval_ratio;
%         results_specs_ratio.exit(j) = exitflag_ratio;
%     catch err
%         disp(tester)        
%     end
% end


%test_eq = (full(Aeq*xstar) -full(beq)) < 1e-3;
j = j + 1;
end
% [e, dL] = simulation_ratio(xstar,p_bar,c(~index_nan_c),w,g,delta,N,10000,0,index_nan_b,index_nan_c, -1);

obj = @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,10000,use_parallel,index_nan_b,index_nan_c, -1); 
tracker = 0
iter = 1
while tracker == 0
    try
        disp(iter)
        knitromatlab( @(x) obj(x),xstar_ratio,Aineq,bineq,Aeq,beq,lb,ub,[],extendedFeatures, options, 'knitro_me_0.opt');%,'knitro.opt'
        tracker = 1
    catch err
        iter = iter + 1
    end
end
% 
% Astar = reshape(xstar(1:N^2),N,N);
% b(index_nan_b)=xstar(N^2+sum(index_nan_b)+1:end);
% c(index_nan_c)=xstar(N^2+sum(index_nan_b)+1:end);
% 
% Astar(Astar<1e-3)=0;
% Gstar = digraph(Astar,names(1:N)); % digraph
% p2=plot(Gstar,'Layout','circle');%'layered','circle','force','subspace','force3','subspace3'
% Gstar.Edges.LWidths =7*Gstar.Edges.Weight/max(Gstar.Edges.Weight);
% p2.LineWidth = Gstar.Edges.LWidths;
% p2.MarkerSize = 30*w/max(w);

%% Use matlab global search 
%'FiniteDifferenceStepSize',0.001,,'StepTolerance',1e-4,'XTolerance',0.1,'ConstraintTolerance',0.001,'OptimalityTolerance',0.01

% options = optimoptions(@fmincon,'Display','iter-detailed','FunctionTolerance',0.01,...
%     'MaxFunctionEvaluations', 10000, 'MaxIterations', 500,...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,...%'FinDiffType','central',...
%     'Algorithm','sqp','ScaleProblem','obj-and-constr');%,'TypicalX',x0+10^(-7));
% problem = createOptimProblem('fmincon','objective',...
%     @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),...
%     'x0',x0,'Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nonlcon(x),'options',options);
% gs = GlobalSearch('Display','iter','FunctionTolerance',0.05);
% 
% [xstar,fval,exitflag,output,solutions] = run(gs,problem);
% % reconstruct A
% Astar = reshape(xstar(1:N^2),N,N);
% % reconstruct b
% b(index_nan_b)=xstar(N^2+sum(index_nan_b)+1:end);
% % reconstruct c
% c(index_nan_c)=xstar(N^2+sum(index_nan_b)+1:end);
% 
% %% Use matlab multistart 
% %,'ScaleProblem','obj-and-constr''StepTolerance',0.001,,'FunctionTolerance',0.1,'ConstraintTolerance',0.001,'OptimalityTolerance',0.001'FiniteDifferenceStepSize',0.001,
% options = optimoptions(@fmincon,'Display','iter-detailed',...
%     'MaxFunctionEvaluations', 100000, 'MaxIterations', 500,...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,...%'FinDiffType','central',...
%     'Algorithm','interior-point',...
%     'FunctionTolerance',0.1,'ConstraintTolerance',0.001,'OptimalityTolerance',0.001);%,'TypicalX',x0+10^(-7)); 
% 
% problem = createOptimProblem('fmincon','objective',...
%     @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),'x0',x0,'Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nonlcon(x),'options',options);
% % problem = createOptimProblem('fmincon','objective',...
% %     @(x) zero_obj(x),'x0',x0,'Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',@(x) nonlcon(x),'options',options);
% 
% ms = MultiStart('Display','iter','UseParallel',false);%,'StartPointsToRun','bounds','FunctionTolerance',0.1,'XTolerance',0.1
% 
% [xstar,fval,exitflag,output,solutions]  = run(ms,problem,20);
% 
% % reconstruct A
% Astar = reshape(xstar(1:N^2),N,N);
% % reconstruct b
% b(index_nan_b)=xstar(N^2+sum(index_nan_b)+1:end);
% % reconstruct c
% c(index_nan_c)=xstar(N^2+sum(index_nan_b)+1:end);
% 
% %%
% sol_num=1;
% Astar = reshape(solutions(sol_num).X(1:N^2),N,N);
% % reconstruct b
% b(index_nan_b)=solutions(sol_num).X(N^2+sum(index_nan_b)+1:end);
% % reconstruct c
% c(index_nan_c)=solutions(sol_num).X(N^2+sum(index_nan_b)+1:end);
% 
% 
% Astar(Astar<1e-3)=0;
% Gstar = digraph(Astar,names(1:N)); % digraph
% p2=plot(Gstar,'Layout','circle');%'layered','circle','force','subspace','force3','subspace3'
% Gstar.Edges.LWidths =7*Gstar.Edges.Weight/max(Gstar.Edges.Weight);
% p2.LineWidth = Gstar.Edges.LWidths;
% p2.MarkerSize = 30*w/max(w);
% 
% p2.NodeColor = repmat([1 0 0],N,1)
% % p2.NodeColor(type'==5,:) = repmat([1 0 0],sum(type==5),1);
% p2.NodeColor(type'==2,:) = repmat([1 1 0],sum(type==2),1);
% p2.EdgeAlpha=0.3
% 
% %% Use matlab pattern search 
% 
% options = optimset('Display','iter','TolX',1e-4,'TolFun',1e-2,'DiffMinChange',0.001);
% [xstar,fval,exitflag,output] = patternsearch(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
% % reconstruct A
% Astar = reshape(xstar(1:N^2),N,N);
% % reconstruct b
% b(index_nan_b)=xstar(N^2+sum(index_nan_b)+1:end);
% % reconstruct c
% c(index_nan_c)=xstar(N^2+sum(index_nan_b)+1:end);
% 
% 
% %% Use matlab fmincon
% %'algorithm','interior-point','HessianApproximation','lbfgs', 'algorithm', 'active-set' 
%  options = optimoptions(@fmincon,'Display','iter-detailed','StepTolerance',0.001,'FunctionTolerance',0.1,'ConstraintTolerance',0.001,'OptimalityTolerance',0.001,...
%     'FiniteDifferenceStepSize',0.001,'MaxFunctionEvaluations', 10000, 'MaxIterations', 500,...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,...%'FinDiffType','central',...
%     'algorithm','sqp' ,'HonorBounds',false);%,'TypicalX',x0+10^(-7));'Algorithm','sqp'
% 
% [xstar,fval,exitflag,output,lambda] = fmincon(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
% % reconstruct A
% Astar = reshape(xstar(1:N^2),N,N);
% % reconstruct b
% b(index_nan_b)=xstar(N^2+sum(index_nan_b)+1:end);
% % reconstruct c
% c(index_nan_c)=xstar(N^2+sum(index_nan_b)+1:end);

%% Close parallel pool, save

% delete(gcp('nocreate'))
save(['sim_results' datestr(now) '_' char(string(num_sim)) 'trials_' char(string(N)) 'nodes'])
diary off