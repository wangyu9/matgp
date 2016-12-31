function [VFs] = process_mesh_sequence(mesh_sequence)

f = size(mesh_sequence,1);
n = size(mesh_sequence{1,1},1);
d = size(mesh_sequence{1,1},2);

VFs = zeros(n,d,f);

for i=1:f
   VFs(:,:,i) = mesh_sequence{i,1}; 
end