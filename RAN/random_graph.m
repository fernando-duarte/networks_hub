function [y,Deg]= random_graph(E,N)
% generates random graph of degree E with N nodes
adj = spalloc(N,N,E);
idx = randperm(N * N, E+N);
idx(ismember(idx, 1:N+1:N*N)) = [];
idx = idx(1:E);
adj(idx) = 1;
adj = min( adj + adj.', 1);
y = full(adj);
Deg=sum(y);

end