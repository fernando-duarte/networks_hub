function [uvector,AD,z] = ui(p_bar,g,x,A,c,N,w)

% solve for clearing vector
y = fsolve(@(p) systemeq(p,p_bar,g,x,A,c),p_bar,optimset('Display','none','TolX',1e-10,'TolFun',1e-10));

%create set of nodes that defaulted
tolerance = 0.001; % value of clearing vector below tolerance is considered equal to zero
z =find(abs(p_bar-y)>tolerance); % set D(x) of defaulting nodes
% z0=setdiff(1:size(A,2),z); % set of non-defaulting nodes
AD = A(z,z);
utemp=(speye(length(z))-(1+g)*AD)\ones(length(z),1);

%create ui
uvector = sparse(1,z,utemp,1,N,length(z));
z2 = find(full(x') > w');
uvector(z2) = max(ones(1, size(z2, 2)), uvector(z2));
% uvector(z)=utemp;


end



