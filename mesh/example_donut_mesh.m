%%
n = 50;
i = (1:n)';
theta = 2*pi/n*i;
VL = [cos(theta),zeros(n,1),sin(theta)];
%%
% scatter(VL(:,1),VL(:,3));
%%
[V,F] = donut_mesh(VL,80,2);
%%
trimesh(F,V(:,1),V(:,2),V(:,3));
%%
writeOBJ('donut.obj',V,F);