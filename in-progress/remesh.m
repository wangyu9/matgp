%%
pdeplot3D(model)

%%
model = createpde(1);
writeSTL('tmp.stl',V,F);
importGeometry(model,'tmp.stl');
%%
tet = generateMesh(model,'GeometricOrder','linear');
%%
pdeplot3D(tet);