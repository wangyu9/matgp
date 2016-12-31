function trisurf_per_element(F,X,Y,Z,AE)

hh = trisurf(F,X,Y,Z);

set(gca, 'CLim', [min(AE), max(AE)]);
set(hh,'FaceColor', 'flat', ...
       'FaceVertexCData', AE, ...
       'CDataMapping', 'scaled');