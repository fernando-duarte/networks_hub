function F = systemeq(p,p_bar,g,x,A,c)

F = p-min(p_bar,max((1+g)*(p*A+c-x)-g*p_bar,0));

end