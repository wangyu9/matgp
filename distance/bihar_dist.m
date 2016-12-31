function [D] = bihar_dist(V,F,b,dim,p)

D_biharmonic = zeros(size(V,1),length(b));
for i=1:size(D_biharmonic,2)
    tmp = biharmonic_distance(V,F,b(i),dim,p);
    D_biharmonic(:,i) = tmp;
end
D = D_biharmonic;