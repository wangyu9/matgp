function [Lz,E] = facet_laplacian_pervertex(V,TF)
  % FACET_LAPLACIAN_PERVERTEX Builds an "edge-based" Laplacian L, which is an #V by #V
  % square matrix which maps scalar functions living at vertices to
  % Laplacian values living at edges. For tets, edges are now facets.
  %
  % [L,E] = facet_laplacian(V,F)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by element-size list of triangle indices
  %   B  #E list bools. true iff unoriented facet occurs exactly once in F
  %     (non-manifold and non-existant edges will be false)
  % Outputs:
  %   L  #V by #V edge-based Laplacian
  %   E  #E by 2 list of edges

[Lez,E] = facet_laplacian(V,TF);
B = is_boundary_facet(E,TF);

% Construct mass matrix
[M,mE] = crouzeix_raviart_massmatrix(V,TF);
% Be sure same edges are being used.
assert(all(E(:)==mE(:)));
% "Kill" boundary edges
Lez(B,:) = 0;

n = size(V,1);
e = size(E,1);


AI = repmat( (1:e)', 1, size(E,2) );
AI = AI(:);
AJ = E(:);
AV = 0.5*ones(size(AI));

A = sparse(AI,AJ,AV,e,n);
%A = bsxfun(@ldivide,sum(A,2),A);

Lz = A'*Lez;