function [c, ceq] = sums(x,N,num_missing)
% passes the non-linear constraints for the optimization
%system 

if any(isnan(x))
    1
end

% reconstruct A
A = reshape(x(1:N^2),N,N);

%  constraints a_ij * a_ji = 0
transpose_A = transpose(A);
ceq_mat = A.*transpose_A;
upper_triang = ~tril(ones(size(ceq_mat)),0);
ceq1 = ceq_mat(upper_triang) ;

c = [];
ceq = ceq1(:);
GC = [];

[ix_row,ix_column]=find(upper_triang);
ix_transpose = fliplr([ix_row ix_column]);

ix =sortrows([sub2ind([N N],ix_transpose(:,1),ix_transpose(:,2)) find(upper_triang)],2);
%GCeq = zeros(N^2,length(ceq));
GCeq = sparse(N^2,length(ceq));

for i=1:length(ceq)
     GCeq(ix(i,:),i) = A(fliplr(ix(i,:)));
end
     
GCeq = [GCeq;sparse(num_missing,N*(N-1)/2)];

end

