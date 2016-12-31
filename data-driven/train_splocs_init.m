%function [D,A] = train_splocs_init(VFs,F,K)


%% VFs (n,3,nf)

%% to remove later
K = 50;%20;
%% input parameter setup
smooth_min_dist = 0.1;
smooth_max_dist = 0.7;
num_iters_max = 10;
sparsity_lambda = 2;

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
%% this ordering is different from splocs paper.
D = zeros(n*dim,K); 
A = zeros(K,nf);
RFs = aVFs;
%% Initializaiton

for k = 1:K
    %%
    % find the vertex explaining the most variance
    magnitude = sum(sum( RFs.^2,3),2);
    [~,idx] = max(magnitude);
    % find the linear component explaining the motion of the vertex
    idx
    % the following SVD is trying to decompose a dim by nf matrix AA(dim,nf)
    % in to AA = UU*SS*VV', or AA' = VV*SS'*UU', so VV(:,1)*SS(1,1) could
    % be interpreted as activation.
    [~,ss,vv] = svd(permute(RFs(idx,:,:),[2 3 1]));% get the 3 by nf parts % equivalent to if remove the first 1-dim: svd(VFs(idx,:,:)');
    Vt = vv';
    % ak is what in the paper as wk. ak should be understood as activations of the new added component for all frames.
    ak0 = ss(1,1)*vv(:,1);% equivalent to: ss(1,1)*Vt(1,:)';
    % ak is nf by 1
    ak_proj = project_weight(ak0);
    ak_proj_negative = project_weight(-ak0);
    
    if(norm(ak_proj)>norm(ak_proj_negative))
        ak = ak_proj;
    else
        ak = ak_proj_negative;
    end
    
    [s,dist] = compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
    s = 1 - s;
    
    % solve for optimal component inside support map
    %  tensor_product( dk , wk) ~= R 
    % dk(n*3,1) is what in the paper as ck
    % r(n*3,nf) is the flatten residual from R(n,3,nf); 
    % wr(n*3,nf) is the weighted residual using s.
    wr = bsxfun(@times,RFs,s);
    wr = reshape(wr,[n*dim,nf]);
    %  dk(n*3,1) * ak(nf,1)' ~= wr(n*3,nf)
    dk = (wr*ak)/(ak'*ak); %ak'*ak is a scale value.
    
    D(:,k) = dk;
    A(k,:) = ak';
            
    RFs  = RFs - reshape(dk*ak',[n,dim,nf]);
    
    % get current reconstruction reuslt. could do it only at the last k
    % this is wrong yet... sploc_init_VFs = reshape(D*A,[n,dim,nf]);
end
D = D / pre_scale_factor;
A_init = A;
D_init = D;
% scale_init = pre_scale_factor;
% shift_init = meanV;

% r(n*3,nf) is the flatten residual from R(n,3,nf); 

%%
if(false)
%% Debug Command for SPLOC Init    
%% to see magnitude.
render_mesh(meanV,F,'ScaleColor',magnitude/max(abs(magnitude)));
%% to see the dist
render_mesh(meanV,F,'ScaleColor',dist/max(dist));
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
render_component(meanV,F,D(:,30),'magnitude',2,'iterations',1);
%%
render_animation(sploc_init_VFs,F);
%%
python_render_component(D,meanV,F);
end