%finds best source and distance to the best source
% if distance is negative, the best source cannot be found (for example, because the propagation was stopped before it reached this point)
% Danil Kirsanov, 09/2007 
% modified by wangyu slightly to make it work: 3/23/2015

function [source_id, distance] = geodesic_best_source(algorithm, destination)

global geodesic_library;

if nargin == 2
    d = geodesic_convert_surface_points({destination});

    tmp = 1;
    [source_id, tmp1, distance] = calllib(geodesic_library, 'distance_and_source', algorithm.id, d, tmp);
    source_id = source_id + 1;
else                                    %return distances and sources for all vertices
    tmp = libpointer('doublePtrPtr');
    tmp1 = libpointer('longPtrPtr');% here is modified by wangyu % original one tmp1 = libpointer('int32PtrPtr');
    
    [num_vertices, d, s] = calllib(geodesic_library, 'distance_and_source_for_all_vertices', algorithm.id, tmp, tmp1);
    
    setdatatype(d, 'doublePtr', num_vertices);
    distance = d.Value;
    setdatatype(s, 'int32Ptr', num_vertices);
    source_id = s.Value + 1;
end