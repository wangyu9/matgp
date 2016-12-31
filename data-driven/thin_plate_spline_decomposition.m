%function [] = thin_plate_spline_decomposition(VFs,HFs)
meanV = mean(VFs,3);
%% Debugging
if(true)
b = farthest_point_sampling(meanV,[],F,[],50,1);
HFs = VFs(b,:,:);
HFs = HFs + 10*(2*rand(size(HFs))-1);
end
%% 
assert(size(HFs,2)==3);
assert(size(VFs,2)==3);
assert(size(VFs,3)==size(HFs,3));
n = size(VFs,1);
m = size(HFs,1);
nf = size(HFs,3);
XFs = reshape(VFs,[n*3,nf]);
kernel = 'biharmonic';%'triharmonic';%'biharmonic';
%% Init
A = reshape(HFs,[m*3,nf]);
meanH = mean(HFs,3);
meanX = mean(XFs,2);
% rrms = norm(bsxfun(@minus,XFs,meanX));
%%
if(true)
%% Mehtod 0: faster
num_iter_coord_decom = 4;
for it = 1:num_iter_coord_decom
    
    sprintf('Iter: %03d',it);
    
    % Step 1: compute W given meanH
    W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH,'kernel',kernel);
    D = blkdiag(W,W,W);
    
    disp(norm(XFs-D*A));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A))
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
        WFs(:,:,j) = ThinPlateSplineCoordinates(VFs(:,:,j),[],3,'C',HFs(:,:,j),'kernel',kernel);
    end
    W = mean(WFs,3);
    D = blkdiag(W,W,W);
    
    disp(norm(XFs-D*A));
    
    % Step 2: update H using current W
    % consider the least square problem XFs ~= D*A
    % (D'*D)*A ~= D'*XFs
    % A ~= (D'*D)\(D'*XFs);
    A = (D'*D)\(D'*XFs);
    HFs = reshape(A,[m,3,nf]);
    meanH = mean(HFs,3);
    
    disp(norm(XFs-D*A))
end
end

%% Debug Commands
if(false)
%%
disp(norm(XFs-D*A)*1000/sqrt(3*nf*n));% E_rms
%%
render_animation(D*A,F);
%%
render_animation(XFs,F);
%%
simple_deform_3D([],meanV,F,meanH,W,'InterpMode','LI');
%%
point_based_skinning(meanV,[],F,W,meanH);
end