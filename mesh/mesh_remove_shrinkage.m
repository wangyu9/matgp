function [V2,F2] = mesh_remove_shrinkage(V,F,Condition,vd)

toRemove = find(~Condition);
toKeep = find(Condition);
%%
n = size(V,1);

r = size(toRemove,1);
k = n - r;

% MAP will map all 
MAP = (1:n)';
MAP(toRemove) = n+1;

MAP2 = -ones(n+1,1);
MAP2(toKeep,1) = (1:k)';
MAP2(n+1,1) = k+1;
%%
F2 = [ MAP(F(:,1)), MAP(F(:,2)), MAP(F(:,3)) ];
%%
FtoRemove = find((F2(:,1)==n+1)+(F2(:,2)==n+1)+(F2(:,3)==n+1)>=2);
%%
F2(FtoRemove,:) = [];
%%
F2 = [ MAP2(F2(:,1)), MAP2(F2(:,2)), MAP2(F2(:,3)) ];

%%
V2 = [V(toKeep,:);vd];


end