function [L,dL_dx]= simulation(x_solve,p_bar,c_non_missing,w,g,delta,N,num_sim,use_parallel,index_nan_b,index_nan_c, min_max)


% reconstruct A
A = reshape(x_solve(1:N^2),N,N);

% reconstruct c
c = sparse(1,N);
c(~index_nan_c)=c_non_missing;
c(index_nan_c)=x_solve(N^2+sum(index_nan_b)+1:end);

% initialize variables
total_losses = zeros(num_sim,1);

num_out = nargout;
if num_out>=2
    %     dL_dA_sim = zeros(num_sim,N^2);
    %     dL_dc_sim = zeros(num_sim,sum(index_nan_c));
    %dL_dx_sim = sparse(num_sim,length(x_solve));
    dL_dx_sim = sparse(num_sim,length(x_solve));
end

% draw random variables
rng('default')

% uniform random variable
%    x = random('Uniform',0,1,[num_sim,N]).*repmat(c,num_sim,1);
u = random('Uniform',0,1,[num_sim,N]);

sim_c = repmat(c,num_sim,1);
sim_w =repmat(w,num_sim,1);
sim_delta=repmat(delta,num_sim,1);
beta_b = log(sim_delta)./log(1-(sim_w./sim_c));
beta_x = betainv(u,1,beta_b);
x =beta_x.*sim_c;
c_dBetaInv_dc=-(((1-u).^(1./beta_b)).*sim_w.*log(1-u).*(1-(sim_w./sim_c)).^(-1) )./(sim_c.*log(sim_delta));
clearvars sim_c sim_w sim_delta beta_b
% c_dBetaInv_dc = zeros(num_sim,N);

% independent beta random variables
% beta_a =1;% ones(1,N);
% if all(beta_a==1)
%     % can solve analytically
%     sim_c = repmat(c,num_sim,1);
%     sim_w =repmat(w,num_sim,1);
%     sim_delta=repmat(delta,num_sim,1);
%    beta_b = log(sim_delta)./log(1-(sim_w./sim_c));
%    beta_x = betarnd(beta_a,beta_b);
%    x =beta_x.*sim_c;
%    dfBeta_dc_fBeta = -sim_w.*(log(1-(sim_w./sim_c))+log(1-beta_x).*log(sim_delta))./(sim_c.*(sim_c-sim_w).*(log(1-(sim_w./sim_c))).^2);
% else
%     beta_b = fsolve(@(b_param) 1-betacdf(w./c,beta_a,b_param)-delta,2*ones(1,N),optimoptions('fsolve','Display','none','TolX',1e-9,'TolFun',1e-9));
% end

%  x = betarnd(1,30,[num_sim,N]).*repmat(c,num_sim,1);

% dependent beta random variables using gaussian copula
% beta_a = ones(1,N); % beta distribution parameter a for each node
% wc_ratio = w./c;
% wc_ratio(wc_ratio>=1) = 0.98;
% [beta_b,beta_b_fval] = fsolve(@(b_param) 1-betacdf(wc_ratio,beta_a,b_param)-delta,ones(1,N),optimoptions('fsolve','Display','iter','TolX',1e-9,'TolFun',1e-9,'Display','none')); % beta distribution parameter b for each node so that delta=prob of default
% if max(abs(beta_b_fval))>1e-3
%     warning('did not find parameters for the beta distribution')
% end
%
% kendall_tau = eye(N);%gallery('randcorr',ones(N,1));%eye(N);%ones(N,N); % pick Kendall's rank correlation
% rho_copula = copulaparam('Gaussian',kendall_tau); % gives linear correlation for gaussian copula
% uniform_x = copularnd('gaussian',rho_copula,num_sim); % create uniformly distributed random variables using gaussian coupla
% beta_x = betainv(uniform_x,repmat(beta_a,num_sim,1),repmat(beta_b,num_sim,1)); % convert so that each marginal has the right beta distribution
% x =beta_x.*repmat(c,num_sim,1); % scale by c

% compute losses for each realization of x
dL_db_sim = sparse(sum(index_nan_b),1);
if use_parallel
    % add derivatives with respect to b and c
    %dL_db_sim = sparse(sum(index_nan_b),1);
    parfor s=1:num_sim
        %         x = random('Uniform',0,1,[1,N]).*c;
        if num_out<2
            total_losses(s)= sum( min(x(s,:),w)+(x(s,:)-w).* ui(p_bar,g,x(s,:),A,c,N));
        else
            %             slice_x = x(s,:);
            %             slice_c_dBetaInv_dc = c_dBetaInv_dc(s,:);
            
            [uvector,AD,z] = ui(p_bar,g,x(s,:),A,c,N);
            total_losses(s)= sum( min(x(s,:),w)+(x(s,:)-w).*uvector );
            
            % create dL/dA
            Id = speye(N); selectD = Id(:,z);
            Id_D = selectD'*Id*selectD;
            uvector_D = uvector*selectD;
            D = length(z);
            left_du_dA = ((Id_D-(1+g)*AD)\Id_D)';
            xw_left_du_dA = repmat(x(s,z)-w(z),D,1).*left_du_dA;
            dL_dA_mat_D = (1+g)* sum(reshape(kron(full(xw_left_du_dA),full(uvector_D)),D,D,D),3);
            dL_dA_mat=selectD*sparse(dL_dA_mat_D)*selectD';
            
            %             temp = false(1,sum(index_nan_c));
            %             temp = zeros(1,N);
            %             for q=1:N
            %                 if x(s,:)<w
            %                     temp(q) = 1;
            %                 end
            %             end
            %                         dL_dc_sim = ((1+g)*uvector(index_nan_c)'+temp).*(slice_x(index_nan_c)./c(index_nan_c))';
            %
            %             dL_dc_sim =(temp(index_nan_c)+uvector(index_nan_c)).*(c_dBetaInv_dc(s,index_nan_c)+x(s,index_nan_c)./c(index_nan_c)  );
            ix_tmp = find(index_nan_c);
            dL_dc_sim =(double(x(s,ix_tmp)<w(ix_tmp))+uvector(ix_tmp)).*(c_dBetaInv_dc(s,ix_tmp)+x(s,ix_tmp)./c(ix_tmp)  );
            
            % put all together
            dL_dx_sim(s,:) = [ dL_dA_mat(:);dL_db_sim;dL_dc_sim(:)];
            %             dL_dx_sim(s,:) = [ dL_dA_mat(:)];
            
        end
        
    end
    %disp('Did One')
else
    
    for s=1:num_sim
        
        if num_out<2
            total_losses(s)= sum( min(x(s,:),w)+(x(s,:)-w).*ui(p_bar,g,x(s,:),A,c,N) );
            if isnan(total_losses(s))
                1
            end
        else
            [uvector,AD,z] = ui(p_bar,g,x(s,:),A,c,N);
            total_losses(s)= sum( min(x(s,:),w)+(x(s,:)-w).*uvector );
            
            % create dL/dA
            Id = speye(N); selectD = Id(:,z);
            Id_D = selectD'*Id*selectD;
            uvector_D = uvector*selectD;
            D = length(z);
            left_du_dA = (1+g)*((Id_D-(1+g)*AD)\Id_D)';
            xw_left_du_dA = repmat(x(s,z)-w(z),D,1).*left_du_dA;
            dL_dA_mat_D =  sum(reshape(kron(full(xw_left_du_dA),full(uvector_D)),D,D,D),3);
            
            dL_dA_mat=selectD*sparse(dL_dA_mat_D)*selectD';
            
            
            ix_tmp = find(index_nan_c);
            dL_dc_sim =(double(x(s,ix_tmp)<w(ix_tmp))+uvector(ix_tmp)).*(c_dBetaInv_dc(s,ix_tmp)+x(s,ix_tmp)./c(ix_tmp)  );
            
            % put all together
            dL_dx_sim(s,:) = [ dL_dA_mat(:);dL_db_sim;dL_dc_sim(:)];
            
            
            %check derivatives
%                                      h =1e-9;
%                                      error_J=nan(N,N);
%                                     for qq=1:N
%                                            for kk=1:N
%                                     J=zeros(N);J(qq,kk)=h;
%                                     [uvector_plus,AD_plus,z_plus] = ui(p_bar,g,x(s,:),A+J,c,N);
%                                     [uvector_minus,AD_minus,z_minus] = ui(p_bar,g,x(s,:),A-J,c,N);
%                                     total_losses_plus= sum( min(x(s,:),w)+(x(s,:)-w).*uvector_plus );
%                                     total_losses_minus= sum( min(x(s,:),w)+(x(s,:)-w).*uvector_minus );
% %                                     [(total_losses_plus-total_losses_minus)/(2*h);   dL_dA_mat(3,3)];
%                                     error_J(qq,kk)=abs( ((total_losses_plus-total_losses_minus)/(2*h)- dL_dA_mat(qq,kk)));
%                                            end
%                                     end
%                                     if abs( ((total_losses_plus-total_losses_minus)/(2*h)- dL_dA_mat(3,3)))  >1e-3
%             
%                                         1
%                                     end
            %
            %                         h =1e-10;
            %                         for p=1:sum(index_nan_c)
            %                             pick_c = p;
            %                             ixc = find(index_nan_c);
            %                             which_c = ixc(pick_c);
            %                             cplus = c; cplus(which_c)=c(which_c)+h;
            %                             cminus = c; cminus(which_c)=c(which_c)-h;
            %                             [uvector_plus,AD_plus,z_plus] = ui(p_bar,g,cplus.*x(s,:)./c,A,cplus,N);
            %                             [uvector_minus,AD_minus,z_minus] = ui(p_bar,g,cminus.*x(s,:)./c,A,cminus,N);
            %                             total_losses_plus= sum( min(cplus.*x(s,:)./c,w)+(cplus.*x(s,:)./c-w).*uvector_plus );
            %                             total_losses_minus= sum( min(cminus.*x(s,:)./c,w)+(cminus.*x(s,:)./c-w).*uvector_minus );
            %
            %                             temp = zeros(sum(index_nan_c),1);
            %                             temp(x(s,index_nan_c)<w(index_nan_c)) = 1; %c(index_nan_c).*
            %                             dL_dc_sim =(temp+(1+g)*uvector(index_nan_c)).*(slice_c_dBetaInv_dc(index_nan_c)+slice_x(index_nan_c)./c(index_nan_c));
            %
            %                             [(total_losses_plus-total_losses_minus)/(2*h);   dL_dc_sim(s,pick_c)]
            %                         end
            
        end
        
    end
end

if any(isnan(total_losses(:)))
    1
end

% compute average losses
L=full(min_max*mean(total_losses,1));

% compute gradient
if num_out>=2
    dL_dx=full(min_max*mean(dL_dx_sim,1));
    
    if max(dL_dx)>1000
        1
    end
end


end
