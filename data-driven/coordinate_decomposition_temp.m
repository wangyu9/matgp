%function [] = thin_plate_spline_decomposition(VFs,HFs,F)
meanV = mean(VFs,3);
%% Debugging
if(true)
b = farthest_point_sampling(meanV,[],F,[],17,1);
HFs = VFs(b,:,:);
end
%% 
assert(size(HFs,2)==3);
assert(size(VFs,2)==3);
assert(size(VFs,3)==size(HFs,3));
n = size(VFs,1);
m = size(HFs,1);
nf = size(HFs,3);
XFs = reshape(VFs,[n*3,nf]);
%% Init
A = reshape(HFs,[m*3,nf]);
meanH = mean(HFs,3);
xVFs = squeeze( VFs(:,1,:) );
yVFs = squeeze( VFs(:,2,:) );
zVFs = squeeze( VFs(:,3,:) );
%%
if(true)
%% Mehtod 0: faster
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
    Wx = (xVFs*xHFs')/(xHFs*xHFs');
    Wy = (yVFs*yHFs')/(yHFs*yHFs');
    Wz = (zVFs*zHFs')/(zHFs*zHFs');
    W = (xVFs*xHFs'+yVFs*yHFs'+zVFs*zHFs')/(xHFs*xHFs'+yHFs*yHFs'+zHFs*zHFs');
    
    W = ( Wx + Wy + Wz )/3;
    %W = ThinPlateSplineCoordinates(meanV,[],3,'C',meanH);
    %D = blkdiag(W,W,W);
    D = blkdiag(Wx,Wy,Wz);
    
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
        WFs(:,:,j) = ThinPlateSplineCoordinates(VFs(:,:,j),[],3,'C',HFs(:,:,j));
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
disp(norm(XFs-D*A)/sqrt(3*n*nf));
%%
render_animation(D*A,F);
%%
render_animation(XFs,F);
%%
simple_deform_3D([],meanV,F,meanH,W,'InterpMode','LI');
end