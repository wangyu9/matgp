%function [D,A] = train_splocs(VFs,meanV,K,D_init,A_init)

% This implements the SIGGRAPH AISA 2013 paper: 
% Sparse Localized Deformation Components

% Link to the paper pdf: http://gvv.mpi-inf.mpg.de/files/splocs.pdf
% A blog post: http://nuit-blanche.blogspot.com/2013/09/sparse-localized-deformation-components.html
% Github source codes: https://github.com/tneumann/splocs
% this implementation is a "translation" from Python to Matlab by wangyu.

%%
K = 50;
[meanV,~,D_init,A_init] = read_splocs('./h5/face_sploc_init.da.h5');
%[~,~,D_init,A_init] = read_splocs('face_splociter=82lambda=1.da.h5');
%%
%render_component(meanV,F,D_init(:,3),'magnitude',2,'iterations',1);
%%
%render_animation( bsxfun(@plus,D_init*A_init,meanV(:)), F);
%% input parameter setup
smooth_min_dist = 0.1/6;
smooth_max_dist = 0.7/6;
num_iters_max = 10;
sparsity_lambda = 1;

rho = 10.0;
num_admm_iterations = 10;

heat_m = 10; % this is the m used for heat time t = m*h^2.

%%
if(false)
    % this is wrong/differernt from the Python codes
    center = mean(meanV);
    aVFs = bsxfun(@minus,VFs,center); % aligned VFs which is (n,dim,f).
else
    % this is similar to PCA such that it delete the meanV from the
    % input sequences. This seems to be important the method to work,
    % at least in the current parameter settings.
    aVFs = bsxfun(@minus,VFs,meanV);
    pre_scale_factor = 1/std(aVFs(:));
    aVFs = aVFs * pre_scale_factor;
end
% render_animation(aVFs,F);
[n,dim,f] = size(aVFs);
aXFs = reshape(aVFs,[n*dim,f]);% aligned XFs, which is (n*dim,f).

h = mean(mean(edge_lengths(meanV,F)));
%[dist,~,~,~,~,dist_pre] = heat_geodesic(meanV,F,1,heat_m*h^2);
[dist,dist_pre] = heat_geodesic_naive(meanV,F,1);

%%
D = D_init * pre_scale_factor;
A = A_init;
assert(size(D,2)==K);
assert(size(A,1)==K);
%% Iteration 0: 
if(true)
    Lambda = zeros(n*dim,K);
    Lambda_reduced = zeros(n,K);
    it = 0;
    % Pre-step 2: update the regularization term will be used in Step 2.
    for k=1:K
        dk = D(:,k);
        % find the one (idx) with biggest displacement in the component
        [~,idx] = max( sum(reshape(dk.*dk,[n,dim]),2) );
        % disp(idx);
        [support_map,dist] = compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
        Lambda_reduced(:,k) = sparsity_lambda * support_map;
        Lambda(:,k) = repmat(sparsity_lambda * support_map,[dim,1]);
    end
    % Step 3:
    RXFs = aXFs - D*A;
    RVFs = reshape(RXFs,[n,dim,nf]);
    sparsity = sum(sum( Lambda.*sqrt(D.^2) ));
    error = sum(sum((RXFs).^2));
    energy = error + sparsity;
    sprintf('Iteration %03d, Energy=%f, Error=%f\n',0,energy,error)
end
%% The Main Global Optimization

% Step 0: prepare auxiliary variables

% Regularization Term:
Lambda = zeros(n*dim,K);
Lambda_reduced = zeros(n,K); % the dim is thrown away as it is same for all.

% Dual Variables:
U = zeros(n*dim,K); % The dimension corresponds to D. 

for it = 1:num_iters_max
    
    sprintf('Iteration %03d\t',it)
    
    % Step 1: fixing D, finding the activation A.
    % the activation C satisfying the condition that max(C(:,k)) == 1, 
    % as the (3) in SPLOC paper.
    
    r = RXFs;
    for k=1:K % for each component
        Dk = D(:,k);
        Dk_norm = Dk'*Dk;
        if(Dk_norm<=1e-8)
           warning(['The component ',num2str(k),' is almost 0, set its activation to 0.\n']);
           A(k,:) = 0;
           continue;
        end
    
        r = r + Dk*A(k,:); % the residual without adding dk
        opt = (r'*Dk)/Dk_norm; %(Dk'*r)/Dk_norm; % opt is nf by 1
        A(k,:) = project_weight(opt);
        r = r - Dk * A(k,:); % update the residual
    end
    
    
    % error reporting
    RXFs = aXFs - D*A;
    RVFs = reshape(RXFs,[n,dim,nf]);
    sparsity = sum(sum( Lambda.*sqrt(D.^2) ));
    error = sum(sum((RXFs).^2));
    energy = error + sparsity;
    sprintf('Iteration %03d: sub 1, Energy=%f, Error=%f\n',0,energy,error)
    
    
    % Pre-step 2: update the regularization term will be used in Step 2.
    for k=1:K
        dk = D(:,k);
        % find the one (idx) with biggest displacement in the component
        [~,idx] = max( sum(reshape(dk.*dk,[n,dim]),2) );
        % disp(idx);
        [support_map,dist] = compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
        Lambda_reduced(:,k) = sparsity_lambda * support_map;
        Lambda(:,k) = repmat(sparsity_lambda * support_map,[dim,1]);
    end
    
    % ADMM (Step 2)
    % argmin_D_Z ||X-D*A||^2 + Omega(Z), s.t. D - Z = 0;
    % 1) D = argmin_D ||X-D*A||^2 + rho*||D-Z+U||^2 
    % 2) Z = argmin_Z Omega(Z) + rho*||D-Z+U||^2
    % 3) U = U + D - Z;
    % note we use rho (as the python codes) instead of rho/2 (as the paper).
    
    % for 1), it is equivalent to update each row of D independently:
    % D(k,:) = argmin_d ||X(k,:)-d'*A||^2 + rho*||d'-Z(k,:)+U(k,:)||^2,
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
    RXFs = aXFs - D*A;
    RVFs = reshape(RXFs,[n,dim,nf]);
    sparsity = sum(sum( Lambda.*sqrt(D.^2) ));
    error = sum(sum((RXFs).^2));
    energy = error + sparsity;
    
    sprintf('Iteration %03d, Energy=%f, Error=%f\n',it,energy,error)
end
%%
D = D / pre_scale_factor;
%%
if(false)
%% Debug Commands for SPLOC Main
%%
figure;
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',support_map);
%% Notice that hte meanV must be added here!!!
disp(1000*norm(XFs-bsxfun(@plus,D*A,meanV(:)),'fro')/sqrt(3*n*nf));
%%
render_component(meanV,F,D(:,50),'magnitude',2,'iterations',1);
%%
python_render_component(meanV,F,D);
%% to see the Labmda distribution.
render_mesh(meanV,F,'ScaleColor',Lambda_reduced(:,1));
%%
render_animation( bsxfun(@plus,D*A,meanV(:)), F);
end