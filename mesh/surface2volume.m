function [V2,F2] = surface2volume(V,F,d)

% given a surface (un-closed mesh), output a volumetric mesh by shifting a
% along the direction vector d
% Then you can do [V3,T3,F3] = tetgen_with_argu(V2,F2,[],'');

% in case the input (V,F) mesh has the normals pointing to z direction
% when using the right hand rule, the normals of output mesh pointing 
% outwards.

assert(size(d,1)==1);
assert(size(d,2)==3);

[V2,F2] = surface2volume_core(V,bsxfun(@plus,V,d),F);