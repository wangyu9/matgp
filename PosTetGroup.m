function PosTetGroup(V,T,k)
%%
[T_cluster] =  T_clustering_via_kmeans(V,V,T,k);

%%
writeDMAT(['T_cluster_V_k=' num2str(k) '.dmat'],T_cluster);