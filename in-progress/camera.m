%%
A = viewmtx(88,78,25);
n = size(V,1);
V4d = [V,ones(n,1)];
V2d = V4d*A';
x2 = V2d(:,1)./V2d(:,4);
y2 = V2d(:,2)./V2d(:,4);
%%
figure;
%plot(x2,y2)
%
trimesh(F,x2,y2);
axis equal;
%%
figure;
trisurf(F,x2,y2,zeros(n,1));