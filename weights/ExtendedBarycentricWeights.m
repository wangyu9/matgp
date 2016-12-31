function [W,V_TC_index] = ExtendedBarycentricWeights(V,C,T_C)
% Wangyu
% C is input control handles
% T_C is #tetrahedra by 4, is the tetrahedra formed by control handles
% ourput: W is the weights matrix
% T_index is which tet on vertex is inside, if no tet, it return -1.
assert(size(V,2)==3);
assert(size(C,2)==3);
if(size(T_C,2)==0)
    error(['The thetrahedra must be given to cal barycentric weights!\n']);
end
assert(size(T_C,2)==4);


% number of vertices
n = size(V,1);
% number of handles
m = size(C,1);
W = zeros(n,m);
V_TC_index = -ones(n,1);

center = [];

for i=1:n
    in_any_tet = false;
    support_handles = [];
    min_dist = 1e9;
    lamada = [];
    for j=1:size(T_C,1)
       % try every tet
       % calculate barycentric weights for each vertex
       support_handles_try = T_C(j,:);
       assert(length(support_handles_try)==4);
       P_try = [C(support_handles_try',:)';[1,1,1,1];];
       lamada_try = inv(P_try)*[V(i,:)';1;];
       if(lamada_try(1)>=-0.000001&&lamada_try(2)>=-0.000001&&lamada_try(3)>=-0.000001&&lamada_try(4)>=-0.000001)
           in_any_tet = true;
           support_handles = support_handles_try;
           P = P_try;
           lamada = lamada_try;
           V_TC_index(i,1) = j;
           break;
       end
       center = mean(C(support_handles_try',:),1);
       dist = (V(i,:)-center)*(V(i,:)-center)';
       if(dist<min_dist)
           min_dist = dist;
           support_handles = support_handles_try;
       end
    end
    if(~in_any_tet)
       P = [C(support_handles',:)';[1,1,1,1];];
       lamada = inv(P)*[V(i,:)';1;];
    end
    W(i,support_handles') = lamada';
end

% make sure that every vertex belongs to one tet
assert(min(min(V_TC_index))>0);



