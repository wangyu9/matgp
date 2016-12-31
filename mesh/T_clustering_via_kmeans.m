function [T_cluster] =  T_clustering_via_kmeans(W,V,TF,k)

dim = size(TF,2)-1;

if(min(min(TF))==0)
    TF = TF+1;warning('Add 1 to TF to make it 1-indexed');
end
assert(min(min(TF))==1);
assert(max(max(TF))==size(V,1));


W_for_TF = zeros(size(TF,1),size(W,2));
for i=1:1:size(TF,1)
    for c=1:1:dim
        W_for_TF(i,:) = W_for_TF(i,:)+W(TF(i,c),:);
    end
end
W_for_TF = W_for_TF/dim;

[T_cluster]=kmeans(W_for_TF,k);

assert(min(T_cluster)==1);

T_cluster = T_cluster - 1;
warning('Mannually make it 0-indexed')
