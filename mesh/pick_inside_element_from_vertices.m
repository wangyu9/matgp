function [FI] = pick_inside_element_from_vertices(TF,VI)

% TF is T or F of a mesh
% PI is the list of indices of the points set
% the return value I is the indices in TF such that all vertices in I,
% which is TF(I,:), belongs to PI.

n = max(max(TF));

Map = zeros(n,1);

Map(VI,:) = 1;

S = zeros(n,1);

for d=1:size(TF,2)
   S = S + Map(TF(:,d),:); 
end

FI = find(S==size(TF,2));