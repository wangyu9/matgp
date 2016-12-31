%function [D,A] = train_splocs(VFs,meanV,K)

% This implements the SIGGRAPH AISA 2013 paper: 
% Sparse Localized Deformation Components

% Link to the paper pdf: http://gvv.mpi-inf.mpg.de/files/splocs.pdf
% A blog post: http://nuit-blanche.blogspot.com/2013/09/sparse-localized-deformation-components.html
% Github source codes: https://github.com/tneumann/splocs
% this implementation is a "translation" from Python to Matlab by wangyu.


%% VFs (n,3,nf)

%% to remove later
K = 50;%20;
%%
smooth_min_dist = 0.1*130;
smooth_max_dist = 0.7*130;
num_iters_max = 10;
sparsity_lambda = 2;

rho = 10.0;
num_admm_iterations = 10;

center = mean(meanV);
aVFs = bsxfun(@minus,VFs,center); % aligned VFs which is (n,dim,f).
% render_animation(aVFs,F);
[n,dim,f] = size(aVFs);
aXFs = reshape(aVFs,[n*dim,f]);% aligned XFs, which is (n*dim,f).

RFs = aVFs;

h = mean(mean(edge_lengths(meanV,F)));
heat_m = 100; % this is the m used for heat time t = m*h^2.
[~,~,~,~,~,dist_pre] = heat_geodesic(meanV,F,1,heat_m*h^2);

% this ordering is different from splocs paper.
D = zeros(n*dim,K); 
A = zeros(K,nf);
%% Initializaiton

for k = 1:K
    %%
    % find the vertex explaining the most variance
    magnitude = sum(sum( RFs.^2,3),2);
    [~,idx] = max(magnitude);
    % find the linear component explaining the motion of the vertex
    idx
    [U,ss,Vt] = svd(permute(RFs(idx,:,:),[2 3 1]));% get the 3 by nf parts % equivalent to if remove the first 1-dim: svd(VFs(idx,:,:)');
    % ak is what in the paper as wk. ak should be understood as activations of the new added component for all frames.
    ak0 = ( ss(1,1) * Vt(1,:) )';
    % ak is nf by 1
    ak_proj = project_weight(ak0);
    ak_proj_negative = project_weight(-ak0);
    
    if(norm(ak_proj)>norm(ak_proj_negative))
        ak = ak_proj;
    else
        ak = ak_proj_negative;
    end
    
    s = 1 - compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
    
    % solve for optimal component inside support map
    %  tensor_product( dk , wk) ~= R 
    % dk(n*3,1) is what in the paper as ck
    % r(n*3,nf) is the flatten residual from R(n,3,nf); 
    r = bsxfun(@times,RFs,s);
    r = reshape(r,[n*dim,nf]);
    %  dk(n*3,1) * ak(nf,1)' ~= r(n*3,nf)
    dk = (r*ak)/(ak'*ak); %ak'*ak is a scale value.
    %dk = reshape(dk,n,dim);
    
    D(:,k) = dk;
    A(k,:) = ak';
    
    %R  = R - reshape(dk*ak',[n,dim,nf]);
    RFs  = RFs - reshape(r,[n,dim,nf]);
    
    % get current reconstruction reuslt. could do it only at the last k
    sploc_init_VFs = reshape(D*A,[n,dim,nf]);
end
A_init = A;
D_init = D;
% r(n*3,nf) is the flatten residual from R(n,3,nf); 
r = bsxfun(@times,RFs,s);
r = reshape(r,[n*dim,nf]);
%% The Main Global Optimization

% Step 0: prepare auxiliary variables

% Regularization Term:
Lambda = zeros(n*dim,K);
Lambda_reduced = zeros(n,K); % the dim is thrown away as it is same for all.

% Dual Variables:
U = zeros(n*dim,K); % The dimension corresponds to D. 

for it = 1:num_iters_max
    
    sprintf('Iteration %03d\t',it)
    
    % Step 1: fixing W, finding the activation C.
    % the activation C satisfying the condition that max(C(:,k)) == 1, 
    % as the (3) in SPLOC paper.
    
    for k=1:K % for each component
        Dk = D(:,k);
        Dk_norm = Dk'*Dk;
        if(Dk_norm<=1e-8)
           warning(['The component ',num2str(k),' is almost 0, set its activation to 0.\n']);
           A(:,k) = 0;
           continue;
        end

        r = r + dk*ak';
        opt = (r'*Dk)/Dk_norm; %(Dk'*r)/Dk_norm; % opt is nf by 1
        A(k,:) = project_weight(opt);
        r = r - Dk * ak';
    end
    
    % Pre-step 2: update the regularization term will be used in Step 2.
    for k=1:K
        dk = D(:,k);
        % find the one (idx) with biggest displacement in the component
        [~,idx] = max( sum(reshape(dk.*dk,[n,dim]),2) );
        support_map = compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
        Lambda_reduced(:,k) = sparsity_lambda * support_map;
        Lambda(:,k) = repmat(sparsity_lambda * support_map,dim,1);
    end
    
    % ADMM (Step 2)
    % argmin_D_Z ||X-D*A||^2 + Omega(Z), s.t. D - Z = 0;
    % 1) D = argmin_D ||X-D*A||^2 + rho*||D-Z+U||^2 
    % 2) Z = argmin_Z Omega(Z) + rho*||D-Z+U||^2
    % 3) U = U + D - Z;
    % note we use rho (as the python codes) instead of rho/2 (as the paper).
    
    % for 1), it is equivalent to update each row of D independently:
    % D(k,:) = argmin_d ||X(k,:)-d'*A||^2 + rho/2*||d'-Z(k,:)+U(k,:)||^2,
    % where d is a K by 1, column vector
    % d = argmin_d d'*(A*A'+rho*I_K)*d ...
    %     + 2*d'*( rho*U(k,:)' - rho*Z(k,:)' -A*X(k,:)' ) + Const.
    %   = argmin_d d'*admm_A*d + 2*d'*admm_b + Const.
    %   = - admm_A\admm_b
    %   where admm_A = A*A' + rho*I_K, 
    %         admm_b = rho*U(k,:)' - rho*Z(k,:)' - A*X(k,:)'
    %  So we have D' = - admm_A\admm_B
    %   where admm_B = rho*U' - rho*Z' - A*X'
     
    % Step 2:
    Z = D; % init dual variables.
    admm_A = A*A' + rho*eye(K);
    admm_B0 = - A * aXFs';
    beta = 1/rho;
    for admm_it=1:num_admm_iterations
       % 1) ADMM Step 1: update D (components)
            % this is equivalent to update each row of D independently
            admm_B = rho*(U-Z)' + admm_B0;
            D = ( - admm_A \ admm_B )';
       % 2) ADMM Step 2: update Z (auxiliary variables)
            % this is equivalent to update each entry of Z independently
            P = D + U;
            Z = max(0,1-beta*Lambda./sqrt(P.*P)).*(P); % the proximal operator
       % 3) ADMM Step 3: update U (dual variables)
            % simple
            U = U + D - Z;
    end
    
    % set updated components to dual Z, 
    % this was also suggested in [Boyd et al.] for optimization of
    % sparsity-inducing norms.
    D = Z;
    
    % Step 3:
    r = aXFs - D*A;
    R = reshape(r,[n,dim,nf]);
    sparsity = sum(sum( Lambda.*sqrt(D.^2) ));
    e = norm(r) + sparsity;
    
    sprintf('Iteration %03d, Energy=%f\n',it,e)
end
%%
if(false)
%% Debug Command for SPLOC Init    
%% to see magnitude.
render_mesh(meanV,F,'ScaleColor',magnitude/max(abs(magnitude)));
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',s);
%% to see the new add component.
render_component(meanV,F,dk,'magnitude',2,'iterations',1);
%% to see the wk over frames.
plot(ak);
render_animation(reshape(dk*ak',[n,dim,nf]),F);
%%
%render_animation(bsxfun(@times,VFs,s),F);
render_animation(reshape(r,[n,dim,nf]),F);
%%
render_animation(R,F);
%%
render_component(meanV,F,D(:,13),'magnitude',2,'iterations',1);
%%
render_animation(sploc_init_VFs,F);
%% Debug Commands for SPLOC Main
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',support_map);
end