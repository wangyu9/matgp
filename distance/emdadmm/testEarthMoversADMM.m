%% Load a mesh

[X, T] = readOff('395.off');
T = double(T); % wangyu

b = 1779;
[D] = earth_mover_distance(X,T,b);
%%
render_mesh3(X,T,'ScaleColor',D);