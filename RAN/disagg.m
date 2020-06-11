%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%   This document constructs constraints for disaggregated nodes    %                                        %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Testing Parameters (Uncomment if you want to test isolated)
% clear
% clc
% clear simulation_ratio


% N        = 3;       % Number of original nodes
N_ins      = 2;       % Num of ins
N_oth      = 1;       % Num of oth
N_reit     = 1;       % Num of reit
 
% insur_dummy = 1;
% reit_dummy  = 0;
% other_dummy = 0;
%  
% index_nan_c = randi([0 1], 1,N_o);
% index_nan_b = randi([0 1], 1,N_o);
%% We need to fix some variables

% check to see which nodes show up
insur_dummy = 1-isempty(find(strcmp(names,'Insurance Aggregate')));
reit_dummy  = 1-isempty(find(strcmp(names,'REIT Aggregate')));
other_dummy = 1-isempty(find(strcmp(names,'Other Aggregate')));

% Number of nodes (if we see the aggregate in the sample)
N_ins       = N_ins  * insur_dummy;      % Num of ins - still need to code in from spreadsheet
N_oth       = N_oth  * other_dummy;      % Num of oth - still need to code in from spreadsheet
N_reit      = N_reit * reit_dummy;       % Num of reit -still need to code in from spreadsheet

% Some counts
N_o         = N;
N_a         = insur_dummy + reit_dummy + other_dummy; % Number of aggregate nodes broken up
N           = N_o + N_ins + N_reit + N_oth - N_a;
N_n         = N_o-N_a;              % Total non-aggregate nodes
N_d         = N-N_n;                % Total dis-aggregated nodes

% Where we want to put ins,oth,reit info
L_ins       = N_o-nnz([N_oth,N_reit]);
L_oth       = N_o-nnz([N_reit]);
L_reit      = N_o;

ind_ins     = N-N_ins-N_oth-N_reit;
ind_oth     = N-N_oth-N_reit;
ind_reit    = N-N_reit;

% Save aggregate data
a_ins       = a(insur_agg);
pbar_ins    = p_bar(insur_agg);
w_ins       = w(insur_agg);
c_ins       = c(insur_agg);
d_ins       = d(insur_agg);
delta_ins   = delta(insur_agg);

a_oth       = a(other_agg);
pbar_oth    = p_bar(other_agg);
w_oth       = w(other_agg);
c_oth       = c(other_agg);
d_oth       = d(other_agg);
delta_oth   = delta(other_agg);

a_reit      = a(reit_agg);
pbar_reit   = p_bar(reit_agg);
w_reit      = w(reit_agg);
c_reit      = c(reit_agg);
d_reit      = d(reit_agg);
delta_reit  = delta(reit_agg);

%remove aggregate data
agg_index = [insur_agg, other_agg, reit_agg];


names(agg_index)    = [];
sectors(agg_index)  = [];
tkr(agg_index)      = [];
a(agg_index)        = [];
p_bar(agg_index)    = [];
c(agg_index)        = [];
b(agg_index)        = [];
d(agg_index)        = [];
f(agg_index)        = [];
w(agg_index)        = [];
delta(agg_index)    = [];
delta  = [delta,repmat(delta_ins,1,N_ins), repmat(delta_oth,1,N_oth), repmat(delta_reit,1,N_reit)];

%revise indecies
index_nan_b         = isnan(b);
index_nan_c         = isnan(c);
index_nan           = index_nan_b | index_nan_c;

%% construct Aeq and beq

% inside assets = sum of inflows 
%   d_i = sum_j[(pbar)_ji]
a1 = [kron(eye(N), ones(1,N))]; 
% revise for original nodes missing c/d
%   a_i = sum_j[(pbar)_ji] + c_i
temp = diag([index_nan_c,zeros(1,N_d)]);
temp(:,all(temp==0))=[];
a1 = [a1,zeros(N,sum(index_nan_b)),temp,zeros(N,3*N_d)]; % b/c/f for aggs irrelevant
temp = [zeros(N_n, N_d); -1*eye(N_d)]; %d is unknown for the new nodes

a1 = [a1,temp];

b1 = [d';zeros(N_d,1)];
b1(index_nan_b) = a(index_nan_b);

%make sure inside assets from non-agg node add up to inside assets
block=1-blkdiag(ones(N_n,N_n), ones(N_ins,N_ins), ones(N_oth, N_oth), ones(N_reit,N_reit));

z1 = zeros(size(block));
z2 = zeros(size(block));
z3 = zeros(size(block));

z1(ind_ins+1:ind_oth,:)     = block(ind_ins+1:ind_oth,:);
z2(ind_oth+1:ind_reit,:)    = block(ind_oth+1:ind_reit,:);
z3(ind_reit+1:end,:)        = block(ind_reit+1:end,:);

a2 = [reshape(z1',1,[]); reshape(z2',1,[]); reshape(z3',1,[])];
a2( all(~a2,2), : ) = []; %delete empty rows
a2 = [a2, zeros(N_a,sum(index_nan_b)+sum(index_nan_c)+4*N_d)];

b2 = [d_ins;d_oth;d_reit];

% inside liabilities = sum of outflows 

a3 = repmat(diag(ones(1,N)),1,N);
% revise for original nodes missing c/d
%   p_bar = sum_i[(pbar)_ji] + f_i
temp = diag([index_nan_b,zeros(1,N_d)]);
temp(:,all(temp==0))=[];
a3 = [a3,temp,zeros(N,sum(index_nan_c)),zeros(N,2*N_d)]; % b/c for aggs irrelevant
temp = [[zeros(N_n, N_d); -1*eye(N_d)],zeros(N,N_d)]; %f is unknown for the new nodes

a3 = [a3,temp];

b3 = [f';zeros(N_d,1)];
b3(index_nan_c) = p_bar(index_nan_c);

% inside liabs and outside liabs from non-agg node add up to total liabs

a4 = [reshape(z1,1,[]); reshape(z2,1,[]); reshape(z3,1,[])];
a4( all(~a4,2), : ) = []; %delete empty rows
temp = [ones(1,N_ins), zeros(1,N_oth+N_reit); zeros(1,N_ins), ones(1,N_oth), zeros(1,N_reit); zeros(1,N_ins+N_oth), ones(1,N_reit)];
temp( all(~temp,2), : ) = []; %delete empty rows

a4 = [a4, zeros(N_a,sum(index_nan_b)+sum(index_nan_c)), temp, zeros(N_a,N_d*3)];

b4 = [pbar_ins;pbar_oth;pbar_reit];

% make sure P_i,i=0
a5 = [full(sparse(1:N,sub2ind([N,N],1:N,1:N),1,N,N^2)), zeros(N,sum(index_nan_b)+sum(index_nan_c)+4*N_d)];

b5 = zeros(N,1);

% make sure outside assets adds up
a6 = [zeros(N_a,N^2 + N_d), zeros(N_a,sum(index_nan_b)+sum(index_nan_c)), temp, zeros(N_a,N_d*2)];

b6 = [c_ins;c_oth;c_reit];

% make sure wealth adds up 
a7 = [zeros(N_a,N^2), zeros(N_a,sum(index_nan_b)+sum(index_nan_c)), -1*temp, temp,-1*temp, temp];

b7 = [w_ins;w_oth;w_reit];

% put it all together
Aeq = [a1;a2;a3;a4;a5;a6;a7];
beq = [b1;b2;b3;b4;b5;b6;b7];

%% Build Aineq and bineq
Aineq1 = [zeros(N_d,N^2+sum(index_nan_b)+sum(index_nan_c)), eye(N_d), -1*eye(N_d), eye(N_d), -1*eye(N_d)];
Aineq2 = [zeros(N_d,N^2+sum(index_nan_b)+sum(index_nan_c)), -1*eye(N_d), zeros(N_d,N_d), -1*eye(N_d), eye(N_d)];
Aineq  = [Aineq1;Aineq2];

bineq  = zeros(N_d*2,1); 

%% Upper and Lower bounds
% elements of the P matrix must be less than the original node's p_bar or
% the aggregate node's p_bar
ub1     = [repmat([p_bar, repmat(pbar_ins,1,N_ins), repmat(pbar_oth,1,N_oth), repmat(pbar_reit,1,N_reit)]',N,1)]; 
% bound missing b's for standard nodes
ub2     = ones(sum(index_nan_b),1).*p_bar(index_nan_b)'- 1e-5;
% bound missing c's for standard nodes
ub3     = ones(sum(index_nan_c),1).*a(index_nan_c)'-1e-5;
% bound b's and c's for aggregate nodes
ub4     = [repmat(pbar_ins,1,N_ins), repmat(pbar_oth,1,N_oth), repmat(pbar_reit,1,N_reit) ,repmat(c_ins,1,N_ins), repmat(c_oth,1,N_oth), repmat(c_reit,1,N_reit)]';
% bound f's and d's for aggregate nodes
ub5     = [repmat(pbar_ins,1,N_ins), repmat(pbar_oth,1,N_oth), repmat(pbar_reit,1,N_reit) ,repmat(d_ins,1,N_ins), repmat(d_oth,1,N_oth), repmat(d_reit,1,N_reit)]';

ub      = [ub1;ub2;ub3;ub4;ub5];
x0helper= [ub1./N_d;ub2./2;ub3./2;ub4./N_d;ub5./N_d];
lb      = zeros(size(ub)); 
% lb    = [zeros(size(ub1));zeros(size(ub2))+1e-5;zeros(size(ub3))+1e-5;zeros(size(ub4))];
