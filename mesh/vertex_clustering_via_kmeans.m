%[V_cluster]=vertex_clustering_via_kmeans(W,k)
W = readDMAT('C:\WorkSpace\Visual Studio 2010\PBS_project\test\arma2\W_bc.dmat');
%%
[V_cluster]=kmeans(W,32);
writeDMAT('V_cluster.dmat',V_cluster);
%% for 3D
T = readDMAT('C:\WorkSpace\Visual Studio 2010\PBS_project\test\arma2\tets.dmat');
[V,F] = readOFF('C:\WorkSpace\Visual Studio 2010\PBS_project\test\arma2\mesh.off');
% [V,T,F] = readMESH('C:\WorkSpace\Visual Studio 2010\PBS_project\test\arma2\mesh.mesh');
TF = T+1;warning('Add 1 to T');
%%
W_for_TF = zeros(size(TF,1),size(W,2));
for i=1:1:size(TF,1)
    for c=1:1:4
        W_for_TF(i,:) = W_for_TF(i,:)+W(TF(i,c),:);
    end
end
W_for_TF = W_for_TF/4;
%%
k = 32;%2*size(W_for_TF,2);
[T_cluster]=kmeans(W_for_TF,k);
%%
if(min(T_cluster)==1)
    T_cluster = T_cluster - 1;warning('Mannually make it 0-indexed')
end
%%
writeDMAT('T_cluster.dmat',T_cluster);
writeDMAT('T_cluster_all0.dmat',zeros(size(T_cluster,1),1));