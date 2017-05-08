function [GL] = graph_laplacian(F)

%%
A = [[F(:,1);F(:,2)],[F(:,2);F(:,3)]];
%%
A = sort(A,2);
%%
A = unique(A,'rows');
%%
G = graph(A(:,1),A(:,2));
GL = laplacian(G);