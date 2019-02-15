function [tH] = total_mean_curvature(V,F)
%%
%V = mat.V;
%F = double(mat.F+1);
%%
L = cotmatrix(V,F);
M = massmatrix(V,F,'full');
H = ((M)\(L*V));
N = per_vertex_normals(V,F);
%%
%tH = sum(sqrt(sum(H.^2,2)));
%tH = sum(sum(H.*N,2));
tH = sum(sum((L*V).*N,2))/sum(sum(M));
%%
%return tH