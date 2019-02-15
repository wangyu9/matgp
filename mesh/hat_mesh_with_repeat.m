function [V,F] = hat_mesh_with_repeat(R,m,n)
% the mesh returned has repeating vertices at north and south pole, and along one longitude line;
% should be removed using meshlab.


% R = 0.5;
% m = 80;
% n = 80;
%
[UV,F] = rectangle_mesh(-1,-1,2/m,2/n,m,n);
F = [F(:,1),F(:,3),F(:,2)];
%%
% render_mesh3(UV,F,'EdgeColor',[0,0,0]);

%%
u = UV(:,1);
v = UV(:,2);
assert(max(u)==1);
assert(max(v)==1);
%%
V = R*[ bsxfun(@times,[cos(pi*u),sin(pi*u)],v+1),1+cos(pi/2*(v+1))];
%%
% render_mesh3(V,F,'EdgeColor',[0,0,0]);