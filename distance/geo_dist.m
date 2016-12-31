function [D] = geo_dist(V,F,b)

assert(size(F,2)==3); % so far supports only 2-D manifolds

% wangyu created the wrapper by modifying from the example 1 from Danil Kirsanov, 09/2007 

global geodesic_library;                
geodesic_library = 'geodesic_release'; %'geodesic_debug';      %"release" is faster and "debug" does additional checks
%rand('state', 0);                         %comment this statement if you want to produce random mesh every time

n = size(V,1);    
m = length(b);


mesh = geodesic_new_mesh(V,F);         %initilize new mesh
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

D = zeros(n,m);              %find distances to all vertices of the mesh (actual pathes are not computed)

for i=1:1:m
    vertex_id = b(i);                             %create sources
    
    source_points = {geodesic_create_surface_point('vertex',vertex_id,V(vertex_id,:))};

    geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)

    [source_id, D(:,i)] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1  

end

geodesic_delete;