%%
n = 20;
len = 5
N = 20;
%
i = (0:n)';
theta = pi*(i*1./n-0.5);
VL1 = [cos(theta)+len, zeros(n+1,1), sin(theta)];
%
VL0 = [len*(0:1:N-1)'/N, zeros(N,1),-ones(N,1)];
VL2 = [len*(N-1:-1:0)'/N, zeros(N,1), ones(N,1)];
%
VL = [VL0;VL1;VL2];
%
%scatter(VL(:,1),VL(:,3));
scatter3(VL(:,1),VL(:,2),VL(:,3));

%%
[V,F] = donut_mesh(VL,120,0);
%%
trimesh(F,V(:,1),V(:,2),V(:,3));
%%
writeOBJ('disk.obj',V,F);