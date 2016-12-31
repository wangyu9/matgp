%function ExternalCallList_TetGroup()

%%
k = 2*size(C,1);
if(dim==3)
    TF = T;
    [T_cluster] =  T_clustering_via_kmeans(W_all,V2,TF,k);
else
    TF = F;
    [T_cluster] =  T_clustering_via_kmeans(W,V,TF,k);
end
%%
writeDMAT('T_cluster.dmat',T_cluster);