%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%   This document converts the results from the disagregated        %
%   network and converts it to constraints for disaggregated nodes  %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Format Agg Data
N_d         = (length(xstar)-N^2-sum(index_nan_b)-sum(index_nan_c))/4;
liab_mat    = reshape(xstar(1:N^2),N,N);
b_aggs      = xstar(end-4*N_d+1:end-3*N_d);
c_aggs      = xstar(end-3*N_d+1:end-2*N_d);
f_aggs      = xstar(end-2*N_d+1:end-N_d);
d_aggs      = xstar(end-N_d+1:end);
p_bar_aggs  = [f_aggs+b_aggs]';
assets_aggs = [d_aggs+c_aggs]';
w_aggs      = assets_aggs-p_bar_aggs;
name_aggs   = [repmat({'ins'},1,N_ins), repmat({'oth'},1,N_oth), repmat({'reit'},1,N_reit)];
sector_aggs = [repmat({'Insurance'},1,N_ins), repmat({'Other'},1,N_oth), repmat({'REIT'},1,N_reit)];

%% Desired Outputs
p_bar       = [p_bar,p_bar_aggs];
w           = [w,w_aggs];
a           = [a,assets_aggs];
b           = [b, b_aggs'];
c           = [c, c_aggs'];
d           = [d, d_aggs'];
f           = [f, f_aggs'];
names       = [names, name_aggs];
tkr         = [tkr, name_aggs];
sectors     = [sectors, sector_aggs];

Astar_ratio_worst = liab_mat./p_bar';

