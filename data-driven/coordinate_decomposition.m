%function [] = thin_plate_spline_decomposition(VFs,HFs,F)
meanV = mean(VFs,3);
%% Debugging
if(true)
    np = 50;
    b = farthest_point_sampling(meanV,[],F,[],np,1);
    HFs = VFs(b,:,:);
    HFs = HFs + 0.3*(2*rand(size(HFs))-1);
elseif(true)
    load b1;
    HFs = VFs(b,:,:);
else
    H2 = readDMAT('meanH.dmat');
    [~,b] = min(pdist2(meanV,H2));
    HFs = VFs(b,:,:);
end
%% 
assert(size(HFs,2)==3);
assert(size(VFs,2)==3);
assert(size(VFs,3)==size(HFs,3));
n = size(VFs,1);
m = size(HFs,1);
dim = 3;
K = m * dim;
nf = size(HFs,3);
XFs = reshape(VFs,[n*3,nf]);
h = mean(mean(edge_lengths(meanV,F)));
%% Init
A = reshape(HFs,[m*3,nf]);
meanH = mean(HFs,3);
xVFs = squeeze( VFs(:,1,:) );
yVFs = squeeze( VFs(:,2,:) );
zVFs = squeeze( VFs(:,3,:) );
%%
if(true)
    
% Step 0: prepare auxiliary variables

% Regularization Term:
%Lambda = zeros(n*dim,K);
Lambda_reduced = zeros(n,m); % the dim is thrown away as it is same for all.

% Dual Variables:
U = zeros(n,m); % The dimension corresponds to D.     

%% Mehtod 0: faster
smooth_min_dist = 0.1;
smooth_max_dist = 0.4;

num_iter_coord_decom = 40;
num_admm_iterations = 10;
rho = 1;
sparsity_lambda = 0.001;
for it = 1:num_iter_coord_decom
    
    sprintf('Iter: %03d',it);
    
    % Pre-step 1: update the regularization term.
    %for k=1:K
        %[support_map,dist] = compute_support_map(idx, meanV, F, heat_m*h^2, dist_pre, smooth_min_dist, smooth_max_dist);
        
        td = pdist2(meanV,meanH);
        td = max(td,smooth_min_dist);
        td = min(td,smooth_max_dist);
        support_map = td;%(td - smooth_min_dist) / (smooth_max_dist - smooth_min_dist);
              
        Lambda_reduced = sparsity_lambda * exp(support_map/(25*h));
        %Lambda(:,k) = repmat(sparsity_lambda * support_map,[dim,1]);
    %end
    
    % Step 1: compute W 
    
    xHFs = squeeze( HFs(:,1,:) );
    yHFs = squeeze( HFs(:,2,:) );
    zHFs = squeeze( HFs(:,3,:) );
        
    if(it==1)    
        W = (xVFs*xHFs'+yVFs*yHFs'+zVFs*zHFs')/(xHFs*xHFs'+yHFs*yHFs'+zHFs*zHFs');
    end
    
    %U = W;
    Z = W; % init dual variables.
    admm_A0 =  (xHFs*xHFs'+yHFs*yHFs'+zHFs*zHFs');
    admm_A = admm_A0 + rho*eye(m);
    admm_B0 = - (xVFs*xHFs'+yVFs*yHFs'+zVFs*zHFs')';
    beta = 1/rho;
    for admm_it=1:num_admm_iterations
       % 1) ADMM Step 1: update D (components)
            % this is equivalent to update each row of D independently
            admm_B = rho*(U-Z)' + admm_B0;
            W = ( - admm_A \ admm_B )';
       % ad hoc
            %W(find(W<0)) = 0; % removing negative W.
       % 2) ADMM Step 2: update Z (auxiliary variables)
            % this is equivalent to update each entry of Z independently
            P = W + U;
            Z = max(0,1-beta*Lambda_reduced./sqrt(P.*P)).*(P); % the proximal operator
            %Z = max(0,1-beta*Lambda_reduced./sqrt(P.*P)).*(P).*(P>=0)...
            %    + max(0,1-beta*8*Lambda_reduced./sqrt(P.*P)).*(P).*(P<0); % the proximal operator with different unbalanced two sides
       % 3) ADMM Step 3: update U (dual variables)
            % simple
            U = U + W - Z;
    end
    
    W = Z;
    %W = ( - admm_A \ admm_B0 )';
    D = blkdiag(W,W,W);
    
    disp(norm(XFs-D*A,'fro'));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A,'fro'));
    
    RXFs = XFs - D*A;
    sparsity = sum(sum( Lambda_reduced.*sqrt(W.^2) ));
    error = sum(sum((RXFs).^2));
    energy = error + sparsity;   
    sprintf('Iteration %03d, Energy=%f, Error=%f\n',it,energy,error)
end
elseif(false)
W = zeros(n,m);
%% Method
num_iter_coord_decom = 4;
for it = 1:num_iter_coord_decom
    
    sprintf('Iter: %03d',it);
    
    td = pdist2(meanV,meanH);
    
    % Step 1: compute W given meanH
    %W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    %W = BilaplacianCoordinates(meanV,F,)
    %Wx * HFs(:,1,:) ~= VFs(:,1,:);

    xHFs = squeeze( HFs(:,1,:) );
    yHFs = squeeze( HFs(:,2,:) );
    zHFs = squeeze( HFs(:,3,:) );
    
    QB = xVFs*xHFs'+yVFs*yHFs'+zVFs*zHFs';
    QA = xHFs*xHFs'+yHFs*yHFs'+zHFs*zHFs';
    
    for j=1:n
        W(j,:) = QB(j,:)/(QA+0.0001*diag(td(j,:).^2));      
    end
    
    D = blkdiag(W,W,W);
    
    %W = ( Wx + Wy + Wz )/3;
    %W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    %D = blkdiag(W,W,W);
    
    
    disp(norm(XFs-D*A,'fro'));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A,'fro'))
end

elseif(true)

%% Method Naive    
num_iter_coord_decom = 4;
for it = 1:num_iter_coord_decom
    
    sprintf('Iter: %03d',it);
    
    % Step 1: compute W given meanH
    %W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    %W = BilaplacianCoordinates(meanV,F,)
    %Wx * HFs(:,1,:) ~= VFs(:,1,:);

    xHFs = squeeze( HFs(:,1,:) );
    yHFs = squeeze( HFs(:,2,:) );
    zHFs = squeeze( HFs(:,3,:) );
    % Wx = (xVFs*xHFs')/(xHFs*xHFs');
    % Wy = (yVFs*yHFs')/(yHFs*yHFs');
    % Wz = (zVFs*zHFs')/(zHFs*zHFs');
    % D = blkdiag(Wx,Wy,Wz);
    
    [~,tps_system] = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    
    W1 = (xVFs*xHFs'+yVFs*yHFs'+zVFs*zHFs')/(xHFs*xHFs'+yHFs*yHFs'+zHFs*zHFs');
    W2 = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    
    W = 1.00 * W1 - 0.00 * W2;
    W(find(W<0)) = 1 * W(find(W<0));
    
    D = blkdiag(W,W,W);
    
    %W = ( Wx + Wy + Wz )/3;
    %W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    %D = blkdiag(W,W,W);
    
    
    disp(norm(XFs-D*A,'fro'));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A,'fro'))
    
    RXFs = XFs - D*A;
    error = sqrt(sum(sum((RXFs).^2)));  
    sprintf('Iteration %03d, Error=%f\n',it,error)
end

else
%% Mehtod 1
num_iter_coord_decom = 4;
WFs = zeros(n,m,nf);
W = mean(WFs,3);
for it = 1:num_iter_coord_decom
    
    sprintf('Iter: %03d',it);
    
    % Step 1: compute W given meanH
    for j = 1:nf
        WFs(:,:,j) = ThinPlateSplineCoordinates(VFs(:,:,j),[],3,'C',HFs(:,:,j));
    end
    W = mean(WFs,3);
    D = blkdiag(W,W,W);
    
    disp(norm(XFs-D*A,'fro'));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A,'fro'))
end
end


%% Debug Commands
if(false)
%%
simple_deform_3D([],meanV,F,meanH,W,'InterpMode','LI','ShowWeightVisualization');
%%
disp(1000*norm(XFs-D*A,'fro')/sqrt(3*n*nf));
%%
render_component(meanV,F,D(:,6+np*0),'magnitude',2,'iterations',1);
%%
render_animation(D*A,F);
%%
render_animation(XFs,F);
%%
compare_animation(XFs,F,D*A,F);
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',support_map(:,1));
%%
simple_deform_3D([],meanV,F,meanH,W,'InterpMode','LI');
%% to see the support distribution.
render_mesh(meanV,F,'ScaleColor',td(:,9));
%%
point_based_skinning(meanV,[],F,W,meanH);
%%
python_render_component(meanV,F,D);
end