function A0 = generateA0(N,row_totals,col_totals)

% temp_A0 = f./p_bar;
%A0 = ones(N)/2;%[0 0 0 0 0 ; 0 0 0 0 0;temp_A0(1:2) 0 temp_A0(4:5); 0 0 0 0 0; 0 0 0 0 0];
% A01= mean((f./p_bar))* random_graph(10,N) ;%triu(ones(N),1)/2;% [ f./p_bar; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
% A02= mean((f./p_bar))*random_graph(10,N) ;%triu(ones(N),1)/2;% [ f./p_bar; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
% A0= triu(A01,1)/2+tril(A02,-1)/2;
% 
% A0 = repmat(((f./p_bar)'/(N-1)),1,N);
% A0(sub2ind([N,N],1:N,1:N))=0;

% A01 = randfixedsum(N,N,min(f./p_bar),0,max(f./p_bar));
% A01(sub2ind([N,N],1:N,1:N))=0;
% A01=triu(A01,1);
% A02 = randfixedsum(N,N,min(d./p_bar),0,max(d./p_bar));
% A02(sub2ind([N,N],1:N,1:N))=0;
% A02=triu(A02,1);
% A0 = (A01+A02')/2;

row_totals = row_totals(:)';
col_totals = col_totals(:);
 
M = N-1;
A0 = zeros(M);
diag_ind = sub2ind([M  M],1:M-1,1:M-1);
for q=1:(M-1)^2;
        if ~ismember(q,diag_ind)
            [ind_row,ind_col]=ind2sub([M-1 M-1],q);
            A0(q) = rand*min(row_totals(ind_col),col_totals(ind_row));
            row_totals(ind_col) = row_totals(ind_col) - A0(q);
            col_totals(ind_row) = col_totals(ind_row) - A0(q);
        end
end


A0 = [A0 col_totals(1:N-1)];
A0 = [A0; [row_totals(1:N-1) 0] ];
% 
% Aeq = [ones(1,N-2) zeros(1,3*(N-2)) 1 0; ... % sum of (N-1)th row
% zeros(1,N-2) ones(1,N-2) zeros(1,2*(N-2)) 0 1;... % sum of (N)th row
% zeros(1,2*(N-2)) ones(1,N-2) zeros(1,N-2) 1 0;... % sum of (N-1)th column
% zeros(1,3*(N-2)) ones(1,N-2) 0 1;... % sum of (N-1)th column
% repmat(eye(N-2),1,2) zeros(N-2,2*(N-2)+2) ;... % sum of rows 1 through N-2
%     zeros(N-2,2*(N-2)) repmat(eye(N-2),1,2) zeros(N-2,2)]; % sum of columns 1 through N-2
% 
% beq = [row_totals(N-1);row_totals(N);col_totals(N-1);col_totals(N);row_totals(1:N-2)';col_totals(1:N-2)];
% 
% x0 = zeros(size(Aeq,2),1);
% lb = zeros(size(x0));
% ub =ones(size(x0));% [repmat(row_totals(N-1),N-2,1);repmat(row_totals(N),N-2,1);repmat(col_totals(N-1),N-2,1);repmat(col_totals(N),N-2,1);min(row_totals(N-1),col_totals(N));min(row_totals(N),col_totals(N-1))];
% 
% options = optimset('Display','iter');
% sol = patternsearch(@(x) zero_obj(x),x0,[],[],Aeq,beq,lb,ub,[],options);
% 
% A0 = [A0 sol(2*M+1:3*M) sol(3*M+1:4*M); sol(1:M)'  0 sol(end-1); sol(M+1:2*M)' sol(end) 0];
% 
% 

end
%%

% m=12;
% A=zeros(600,m);
% cind=randperm(size(A,1));
% 
% for n=1:numel(cind)
%         ind=find(sum(A)<53);
%         try
%             A(cind(n),ind(randperm(numel(ind),4)))=1;
%         catch err
%         end
% end

%%
% 
% function [M]=matsrand(n,c)
% 
%     MM=0;   %arbitrary starting value
%     all_nums=1:n;
% 
%     while MM ~=n*c
% 
%         M = sparse([],[],[],n,n,n*c);
%         ctin = zeros(1,n);
% 
%         for ii=1:n
%             noconnect=ctin>=c;
%             noconnect(ii)=true;
% 
%             rem_nums = all_nums(~noconnect); % remaining numbers
%             rp=randperm(n-sum(noconnect));
%             rp = rem_nums(rp); % remaining numbers, hussled
% 
%             if numel(rp)<c
%                 break
%             else
%                 r=rp(1:c);
%             end
%             M(ii,r)=1;
%             ctin(r)=ctin(r)+1;
%         end
%         MM=sum(ctin);
%     end
% end