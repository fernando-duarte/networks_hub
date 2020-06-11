function [cineq, ceq] = nonlincons_system(y,p,g,x,A,c)

cineq = [];
ceq = systemeq(y,p,g,x,A,c);

end 