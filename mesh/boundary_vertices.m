function [B,J] = boundary_vertices(TF)
  % BOUNDARY_VERTICS Determine boundary verticess of triangle/tetrahedra
  % stored in TF
  %
  % F = boundary_vertices(T)
  %
  % Input:
  %   TF  tetrahedron index list, m by 4, where m is the number of triangle/tetrahedra
  % Output:
  %   B  list of boundary vertices, n by 1, where n is the number of boundary vertices
  %   J  list of indices into TF, n by 1
  % wangyu

  if(size(TF,2)==4)
    % 3D
    [B,J] = boundary_faces(TF);
    J = [J,J,J];
  else
    % 2D
    [B,J] = boundary_edges(TF);
    J = [J,J];
  end

  B = B(:);
  J = J(:);
  
  [B,IA] = unique(B,'rows');
  J = J(IA);
  
end