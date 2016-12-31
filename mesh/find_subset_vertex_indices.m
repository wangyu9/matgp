function [J] = find_subset_vertex_indices(V,Vsub)

[CC,IA,IC] = intersect(V,Vsub,'rows');
assert(size(IA,1)==size(Vsub,1));
% We have CC = V(IA,:) = V_input(IC,:) = 
%   sparse(1:size(IA,1),IA,1,size(IA,1),size(V,1))*V = sparse(1:size(IC,1),IC,1)*V_input
% so V_input =
% sparse(1:size(IC,1),IC,1)\sparse(1:size(IA,1),IA,1,size(IA,1),size(V,1)) * V;
[~,J] = find(sparse(1:size(IC,1),IC,1,size(IC,1),size(V,1))\sparse(1:size(IA,1),IA,1,size(IA,1),size(V,1))==1);
% then we have V_input = V(J,:)
assert( max(max(Vsub-V(J,:)))<0.0001 );