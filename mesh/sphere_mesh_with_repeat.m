function [V,F] = sphere_mesh_with_repeat(R,m,n)
% the mesh returned has repeating vertices at north and south pole, and along one longitude line;
% should be removed using meshlab.

% R = 1;
% m = 20;
% n = 20;
%%
[UV,F] = rectangle_mesh(-1,-1,2/m,2/n,m,n);

%%
render_mesh3(UV,F,'EdgeColor',[0,0,0]);

%%
u = UV(:,1);
v = UV(:,2);
%%
assert(max(u)==1);
assert(max(v)==1);

%%
V = R*[cos(pi*u).*cos(pi/2*v),sin(pi*u).*cos(pi/2*v),sin(pi/2*v)];
%%
render_mesh3(V,F,'EdgeColor',[0,0,0]);