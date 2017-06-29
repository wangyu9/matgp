function h = showDescriptor(X, T, f)
% Displays per-face or per-vertex function f on mesh (X,T)

h = figure;
nv = size(X,1);
figure(h);
if size(f,1) == nv
    patch('vertices',X,'Faces',T,'FaceColor','interp','CData',double(f),'edgecolor','none'); 
else
    patch('vertices',X,'Faces',T,'FaceColor','flat','CData',double(f),'edgecolor','none'); 
end

axis equal;
colorbar;
cameratoolbar;
axis off;