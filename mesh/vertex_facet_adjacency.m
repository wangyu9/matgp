function [VT] = vertex_facet_adjacency(F)
  % VERTEX_FACET_ADJACENCY Build a vertex-facet adjacency matrix for a
  % triangle mesh or tet mesh.
  %
  % VT = vertex_triangle_adjacency(F)
  %
  % Input:
  %   F   #F x 3|4   list of facet indices
  % Output:
  %   VT  #F x #V  sparse matrix, VT(i,j) is 1 iff vertex i is in facet j

  i = (1:size(F,1))';
  j = F;
  VT = sparse(repmat(i,[1,size(F,2)]),j,1);

end
