function [D] = heat_kl_dist(V,F,t,b1)


n = size(V,1);

L = cotmatrix(V,F);
%[L,~] = facet_laplacian_pervertex(V,F);
%[~,~,L] = laplacian_and_mass(V,F);

expLt = expm(L*t);

%
KL = @(p,q) p.*log(p./q)+q.*log(q./p) ;
%
%b1 = 1500;
D = bsxfun(KL,expLt(:,b1)',expLt);
D = sum(D,2);

%%
%trimesh(F, V(:,1),V(:,2), D);
%%
%trimesh(F, V(:,1),V(:,2), sqrt(D));