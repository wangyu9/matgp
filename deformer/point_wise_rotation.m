function [U] = point_wise_rotation(V,T)
  %  Compute point wise rotation of vertices V, using
  % rotations at some control points T, propogated to the mesh using
  % weights W.
  %
  % [U] = lbs(V,T,W)
  % 
  % Inputs:
  %  V  list of vertex positions
  %  T  list of transformations for each controls point, for 2D:
  %    2 by 2 by #vertices, for 3D: 3 x 3 by # vertices
  % Output:
  %  U  list of new vertex positions
  %
  
AV = zeros([size(V',1),1,size(V',2)]);
AV(:,1,:) = V';
U = sum( T .* permute( repmat(AV,[1,3,1]), [2,1,3]), 2 );
U = squeeze(U);
U = U';