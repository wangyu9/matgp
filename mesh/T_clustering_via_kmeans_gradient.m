function [T_cluster] =  T_clustering_via_kmeans_gradient(W,V,TF,k)

dim = size(TF,2)-1;

if(min(min(TF))==0)
    TF = TF+1;warning('Add 1 to TF to make it 1-indexed');
end
assert(min(min(TF))==1);
assert(max(max(TF))==size(V,1));

assert(dim==3);
% only 3d is implemented yet

G = zeros(size(TF,1),size(W,2)*dim);
P = zeros(dim,dim);
B = zeros(dim,size(W,2));
for i=1:1:size(TF,1)
    for c=1:1:dim
        P(c,:) = V(TF(i,c+1),:) - V(TF(i,1),:);
        B(c,:) = W(TF(i,c+1),:) - W(TF(i,1),:);
        %W_for_TF(i,:) = W_for_TF(i,:)+W(TF(i,c),:);
    end
    Gi = P\B;
    GiT = Gi';
    G(i,:) = (GiT(:))';
end

[T_cluster]=kmeans(G,k);

assert(min(T_cluster)==1);

T_cluster = T_cluster - 1;
warning('Mannually make it 0-indexed')
