function [D_ring] = cage_ring_dist(V,C,T_C,V_TC_index)
% by wangyu. 
% V is all vertices
% C is all cage vertices
% T_C is all cage tets 
% V_TC_index is telling each vertex belongs to which tets  

% ourput: D_ring is #vertices by #C

addpath ./gaimc

allE = [];

% add faces edges
% allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
% notice each edge is added twice, but it does not matter later

% add tet edges
if(~isempty(T))
    tetE = [T(:,[1 2]);T(:,[2 3]);T(:,[1 3]);T(:,[1 4]);T(:,[2 4]);T(:,[3 4])];
    allE = [allE;tetE];
    % notice each edge is added twice, but it does not matter later
end
% Map duplicate edges to first instance
%[E,~,EMAP] = unique(sort(allE,2),'rows');

n = size(V,1);
m = size(C,1);

E = allE;

%Dist = dist(V');% do not use this, it cal too much

A11 = sparse(n,n);%this should be zero 
A12 = sparse(n,m);
A21 = sparse(m,n);
for k=1:size(E,1)
    i=E(k,1);
    j=E(k,2);
%     v1 = F(i,:);
%     v2 = F(j,:);
%     dist = sqrt();
    assert(i<=n);
    assert(j<=m);
    A12(i,j) = 1;%sqrt(sum((V(i,:) - V(j,:)).^2));
    A21(j,i) = 1;%sqrt(sum((V(i,:) - V(j,:)).^2));
end

% if one vertext is inside the tet, then the distance is set to 1.
A((1:n)'+n*(V_TC_index-1)) = 1;

%A = spfun(@(x) x-min(min(A))+1,A); % remove the negative edges

As = convert_sparse(A);
%
% Now, we'll run Dijkstra's algorithm for every vertex and save the result
% On my 2GHz laptop, this takes 0.000485 seconds.
isequal(n,size(A,1));
D_ring = zeros(n,length(b));

tic
for i=1:length(b)
    D(:,i) = dijkstra(As,b(i));
end
toc
