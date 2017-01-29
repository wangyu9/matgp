%function [] = heat_kernel_signature()
error(['Not yet to publish!']);
%%
p = 600;
t = 1e-4;
TF = F;
%% removing unreferenced vertices here.
M = massmatrix(V,TF,'barycentric');
b = find(sum(M,2) == 0);
V(b,:) = [];

%%
% get cotangent matrix
L = cotmatrix(V,TF);
% This should be better, but yaron seemed to use barycentric
%M = massmatrix(V,F);
M = massmatrix(V,TF,'barycentric');

%%
[EV,ED] = eigs(-L,M,p+1,'sm');

%EV = EV(:, 2:end);
%ED = ED(2:end, 2:end);

%%
B = sum( EV.^2*exp(-t*abs(ED)), 2);
%B = B(:,1);
%%
trisurf(F,V(:,1),V(:,2),B,B);
%% Debug 
if(false)
%%
trisurf(F,V(:,1),V(:,2),B*100,B*100);
%%
trisurf(F,V(:,1),V(:,2),L*B,B*100);
%%
render_mesh(V,F,'ScaleColor',B./max(B)); 
%%
[D] = biharmonic_distance(V,F,100,100,2);
trisurf(F,V(:,1),V(:,2),D,D);
end