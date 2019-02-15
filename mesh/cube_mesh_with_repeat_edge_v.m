function [V,F] = cube_mesh_with_repeat_edge_v(n)
% note the mesh generated has repeating vertices at edges and corner, so a
% removal of duplicated vertices using e.g. meshlab is necessary.
%%

[V0,F0] = symmetric_rectangle_mesh(-1,-1,1/n,1/n,n,n);
V0 = [V0,ones([size(V0,1),1])];
%%
render_mesh3(V0,F0,'EdgeColor',[0,0,0]);
%%
V = V0;
F = F0;
[V,~,F] = merge_mesh(V,[],F,V0*axisangle2matrix([1 0 0],1*pi/2),[],F0);
[V,~,F] = merge_mesh(V,[],F,V0*axisangle2matrix([1 0 0],2*pi/2),[],F0);
[V,~,F] = merge_mesh(V,[],F,V0*axisangle2matrix([1 0 0],3*pi/2),[],F0);
[V,~,F] = merge_mesh(V,[],F,V0*axisangle2matrix([0 1 0], pi/2),[],F0);
[V,~,F] = merge_mesh(V,[],F,V0*axisangle2matrix([0 1 0],-pi/2),[],F0);
