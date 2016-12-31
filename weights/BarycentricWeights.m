function [W] = BarycentricWeights(support_region,V,bc,b)
% Wangyu
% support_region: #vertices by #handles, value is 1 if the vertex is support by the correspond handles
% Note the column of support_region should be in handle INPUT order 
% each vertex must be supported by exact 4 handles
if(length(find(sum(support_region,2)~=4))~=0)
    error(['Error: Some vertex is not supported by exact 4 handles!\n']);
end

[row_b, col_b] = find(bc==1);
bb= b(row_b);% bb is list of handles in INPUT order. b is not in input order.

% number of vertices
n = size(V,1);
% number of handles
m = size(bc,2);
W = support_region*1.0;
for i=1:n
   % calculate barycentic weights for each vertex
   [~,support_handles] = find(support_region(i,:)==1);
   assert(length(support_handles)==4);
   P = [V(bb(support_handles'),:)';[1,1,1,1];];
   lamada = inv(P)*[V(i,:)';1;];
   %lamada = (lamada<100).*lamada.*(lamada>-100) + (lamada>100)*1 + (lamada<-100)*(-100);% threshold too large weight
   W(i,support_handles') = lamada';
end




