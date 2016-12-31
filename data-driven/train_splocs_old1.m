%function [] = train_splocs(VFs,meanV,K)

% This implements the SIGGRAPH AISA 2013 paper: 
% Sparse Localized Deformation Components

% Link to the paper pdf: http://gvv.mpi-inf.mpg.de/files/splocs.pdf
% A blog post: http://nuit-blanche.blogspot.com/2013/09/sparse-localized-deformation-components.html
% Github source codes: https://github.com/tneumann/splocs



%% VFs (n,3,nf)

%% to remove later
K = 50;%20;
%%
smooth_min_dist = 0.1*140;
smooth_max_dist = 0.7*140;
num_iters_max = 10;
sparsity_lambda = 2;

rho = 10.0;
num_admm_iterations = 10;

center = mean(meanV);
aVFs = bsxfun(@minus,VFs,center);
% render_animation(aVFs,F);
[n,dim,f] = size(aVFs);

R = aVFs;

h = mean(mean(edge_lengths(meanV,F)));
heat_m = 100; % this is the m used for heat time t = m*h^2.
[~,~,~,~,~,dist_pre] = heat_geodesic(meanV,F,1,heat_m*h^2);

% this ordering is different from splocs paper.
W = zeros(n*dim,K); 
C = zeros(K,nf);
%% Initializaiton

for k = 1:K
    %%
    % find the vertex explaining the most variance
    magnitude = sum(sum( R.^2,3),2);
    [~,idx] = max(magnitude);
    % find the linear component explaining the motion of the vertex
    idx
    [U,ss,Vt] = svd(permute(R(idx,:,:),[2 3 1]));% get the 3 by nf parts % equivalent to if remove the first 1-dim: svd(VFs(idx,:,:)');
    % ak is what in the paper as wk. ak should be understood as activations of the new added component for all frames.
    ak0 = ss(1,1) * Vt(1,:);
    % ak is 1 by nf
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
    r = bsxfun(@times,R,s);
    r = reshape(r,[n*dim,nf]);
    %  dk(n*3,1) * ak(1,nf) ~= r(n*3,nf)
    dk = (r*ak')/(ak*ak'); %ak*ak' is a scale value.
    %dk = reshape(dk,n,dim);
    
    W(:,k) = dk;
    C(k,:) = ak;
    
    %R  = R - reshape(dk*ak,[n,dim,nf]);
    R  = R - reshape(r,[n,dim,nf]);
    
    % get current reconstruction reuslt. could do it only at the last k
    sploc_init_VFs = reshape(W*C,[n,dim,nf]);
end
C_init = C;
W_init = W;
% r(n*3,nf) is the flatten residual from R(n,3,nf); 
r = bsxfun(@times,R,s);
r = reshape(r,[n*dim,nf]);
%% The Main Global Optimization
% prepare auxiliary variables
Lambda = zeros(K,n);
for it = 1:num_iter_max
    % Step 1: fixing W, finding the activation C.
    % the activation C satisfying the condition that max(C(:,k)) == 1, 
    % as the (3) in SPLOC paper.
    
    for k=1:K % for each component
        Dk = W(:,k);
        Dk_norm = Dk'*Dk;
        if(Dk_norm<=1e-8)
           warning(['The component ',num2str(k),' is almost 0, set its activation to 0.\n']);
           C(:,k) = 0;
           continue;
        end

        r = r + dk*ak;
        opt = (Dk'*r)/Dk_norm; % opt is 1 by nf
        C(:,k) = project_weight(opt');
        r = r - Dk * ak;
    end
    
end
%%
if(false)
%% to see magnitude.
render_mesh(meanV,F,'ScaleColor',magnitude/max(abs(magnitude)));
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',s);
%% to see the new add component.
render_component(meanV,F,ck,'magnitude',2,'iterations',1);
%% to see the wk over frames.
plot(wk);
render_animation(reshape(ck*wk,[n,dim,nf]),F);
%%
%render_animation(bsxfun(@times,VFs,s),F);
render_animation(reshape(r,[n,dim,nf]),F);
%%
render_animation(R,F);
%%
render_component(meanV,F,W(:,8),'magnitude',2,'iterations',1);
%%
render_animation(sploc_init_VFs,F);
end