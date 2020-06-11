%% 
clear
clc

%% 
use_parallel =0; 
pool_size = 60;
if use_parallel
    p = gcp('nocreate'); 
    if isempty(p)
         % if pool not open, create one
         mypool=parpool(pool_size,'AttachedFiles',{'ui.m','systemeq.m','sums.m'});
    end
end

% Number of trials for simulation
num_sim = 1;
%% set up network
y=10;
p_bar = [55+y 55+y 140 55+y 55+y]; % tot liab
A0 = [0 y/p_bar(1) 0 0 0; 0 0 0 0 y/p_bar(1); 1/14 1/14 0 1/14 1/14; y/p_bar(4) 0 0 0 0; 0 0 0 y/p_bar(5) 0];
A0 = A0(:);

c = [50 50 150 50 50]; % outside assets
w = [5 5 10 5 5]; % net worth
g = 0.5;
b = [55 55 100 55 55]; % outside liabilities
a = w+p_bar; % total assets
d=  a-c;% inside assets
f = p_bar-b;% inside liabilities
N = length(c); % number of nodes
delta = 0.2*ones(1,N); % probabilities of default

%% set up missing b and c
index_nan_b = false(size(b)); index_nan_b(4)=true;
index_nan_c = false(size(c)); index_nan_c(4)=true;index_nan_c(5)=true;

index_nan  = index_nan_b | index_nan_c;

%% set up x0 and dx

h=1e-10;
x0 =A0 ;

x0 = [x0;b(index_nan_b)';c(index_nan_c)'];
b(index_nan_b)=NaN;
c(index_nan_c)=NaN;
gr_dx=nan(length(x0),1);

for q=1:length(x0)
variable = q; % position of variable whose derivative we want to check

dx = zeros(size(x0));dx(variable)=h;
xplus = x0+dx;
xminus = x0-dx;



%% compute gradient using function or finite differences

Lplus =simulation(xplus,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c);
Lminus = simulation(xminus,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c);

gr_dx(q) = (Lplus-Lminus)/(2*h);
end

[Los,gr]= simulation(x0,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c);

both = [gr',gr_dx];
max(abs(gr'-gr_dx))

%% check with optmizer
%options=optimoptions('fmincon','Display','iter','MaxIter',1,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',true,'FinDiffType','central');

% non-linear constraints
 %nonlcon = @(x) sums(x,N,sum(index_nan_b)+sum(index_nan_c));
fmincon(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,[],[],[],[],[],[],nonlcon,options);
fmincon(@(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,[],[],[],[],[],[],[],options);

nonlcon = @(x) sums(x,N,sum(index_nan_b)+sum(index_nan_c));
options = optimset('FinDiffType','central', 'GradObj','on','GradConstr','off','DerivativeCheck','on','MaxIter',0);
[xstar,fval,exitflag,output,lambda,grad] =knitromatlab( @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,[],[],[],[],[],[],nonlcon,[],options);



lb_A=sparse(N^2,1);   ub_A = repmat((f./p_bar)',N,1); ub_A(isnan(ub_A))=1;
lb_b=sparse(sum(index_nan_b),1)+1e-5;   ub_b = ones(sum(index_nan_b),1).*p_bar(index_nan_b)'- 1e-5;
lb_c=ones(sum(index_nan_c),1).*w(index_nan_c)'+1e-5;   ub_c = ones(sum(index_nan_c),1).*a(index_nan_c)'-1e-5;

lb = [lb_A;lb_b;lb_c];
ub = [ub_A;ub_b;ub_c];
[xstar,fval,exitflag,output,lambda,grad] =knitromatlab( @(x) simulation(x,p_bar,c(~index_nan_c),w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c),x0,[],[],[],[],lb,ub,[],[],options,'knitro.opt');








% delete(gcp('nocreate'))

