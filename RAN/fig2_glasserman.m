
clear

 
%Network parameters as in all other files
y=10;
p_bar = [55+y 55+y 140 55+y 55+y]; % tot liab
A0 = [0 y/p_bar(1) 0 0 0; 0 0 0 0 y/p_bar(1); 1/14 1/14 0 1/14 1/14; y/p_bar(4) 0 0 0 0; 0 0 0 y/p_bar(5) 0];
A0 = A0(:);

c = [50 50 150 50 50]; % outside assets
w = [5 5 10 5 5]; % net worth
g = 0.05;

% Network primitives
b = [55 55 100 55 55]; % outside liabilities
a = w+p_bar; % total assets
d=  a-c;% inside assets
f = p_bar-b;% inside liabilities
N = length(c); % number of nodes

pij=A0.*repmat(p_bar',1,N)

[sum(A0,2).*p_bar' p_bar'-b']
[sum(A0.*repmat(p_bar',1,N),1); a-c]
[w-c+p_bar; sum(A0.*repmat(p_bar',1,N),1)]
 
shock = [0 0 0 0 0];
p_initial = fsolve(@(p) systemeq(p,p_bar,g,shock,A0,c),p_bar-1);

%%
% linear inequality constraints 
Aineq = [];bineq=[];%Aineq=repmat(eye(N),1,N); bineq = (f./p_bar)'; % sum of columns of A cannot exceed ratio of inside to total liabilities

% linear equality constraints 
Aeq1=kron(eye(N),p_bar); beq1=d';  % for each i, sum( p_bar * ith column of A) = inside assets
Aeq2=repmat(diag(p_bar),1,N); beq2 =f';  %for each i, p_bar_i*sum( ith row of A) = inside liabilities
Aeq3=full(sparse(1:N,sub2ind([N,N],1:N,1:N),1,N,N^2));beq3=zeros(N,1); % diagonal elements of A must be zero

Aeq = [Aeq1;Aeq2;Aeq3];
beq = [beq1;beq2;beq3]; 

% upper and lower bounds: elements of A between 0 and ratio of inside to total liabilities
lb = zeros(size(A0))-0.0001;
ub = repmat((f./p_bar)',N,1)+0.0001;%ones(size(A0));

% non-linear constraints
nonlcon = @(A) sums(A,N);