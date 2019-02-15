R = 1;
m = 20;
n = 20;

[UV,F] = rectangle_mesh(-1,-1,2/m,2/n,m,n);

%
u = UV(:,1);
v = UV(:,2);
assert(max(u)==1);
assert(max(v)==1);
%
V = R*[ bsxfun(@times,[cos(pi*u),sin(pi*u)],cos(pi/2*v)+0.5*(sin(pi/2*v)+1)),sin(pi/2*v)];
%%
render_mesh3(V,F,'EdgeColor',[0,0,0]);
%%
[V2,F2] = meshlab_mlx(V,F,'remove_duplicated_vertex');
%%
[V2,F2] = meshlab_mlx(V,F,'merge_close_vertices');
%%
render_mesh3(V2,F2,'EdgeColor',[0,0,0]);