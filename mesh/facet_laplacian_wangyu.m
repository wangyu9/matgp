function [L,E] = facet_laplacian_wangyu(V,F)
  % FACET_LAPLACIAN Builds an "edge-based" Laplacian L, which is an #E by #V
  % rectangular matrix which maps scalar functions living at vertices to
  % Laplacian values living at edges. For tets, edges are now facets.
  %
  % [L,E] = facet_laplacian(V,F)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by element-size list of triangle indices
  % Outputs:
  %   L  #E by #V edge-based Laplacian
  %   E  #E by 2 list of edges
  %
  % Examples:
  %   % load some non-convex shape (V,F)
  %   [L,E] = facet_laplacian(V,F);
  %   % find boundary edges
  %   B = is_boundary_edge(E,F);
  %   % Construct mass matrix
  %   [M,mE] = crouzeix_raviart_massmatrix(V,F);
  %   % Be sure same edges are being used.
  %   assert(all(E(:)==mE(:)));
  %   % "Kill" boundary edges
  %   L(B,:) = 0;
  %   % Linear functions are now in the spectrum
  %   [~,~,U] = svd(full(L));
  %   tsurf(F,[V(:,1:2) U(:,end-1)])
  %
  %   % mesh in (V,F)
  %   [Le] = facet_laplacian(V,F);
  %   [Lcr,E,EMAP] = crouzeix_raviart_cotmatrix(V,F);
  %   V2E = sparse(E(:),repmat(1:size(E,1),1,2)',1,size(V,1),size(E,1))';
  %   V2E = bsxfun(@rdivide,V2E,sum(V2E,2));
  %   LcrV2E = Lcr*V2E;
  %   max(max(abs(LcrV2E--2*Le)))
  %
  %
  % See also: crouzeix_raviart_massmatrix, is_boundary_edge
  %

  %% check for non-manifold edges
  %S = statistics(V,F,'Fast',true);
  %if S.num_nonmanifold_edges > 0
  %  error(sprintf('There are %d non-manifold edges',S.num_nonmanifold_edges));
  %end

  switch size(F,2)
  case 3
    % number of vertices
    n = size(V,1);
    % number of faces
    m = size(F,1);
    % Compute cotangents
    C = cotangent(V,F);
    % Map each face-edge to a unique edge
    F2E = reshape(1:3*m,m,3);
    % Assemble entries
    LI = [ F2E(:,[1 1 1 2 2 2 3 3 3]) ];
    LJ = [   F(:,[1 2 3 2 3 1 3 1 2]) ];
    % allE 1st    2 2 2 
    % allE 2nd    3 3 3
    LV = [ ...
      C(:,2)+C(:,3), -C(:,3), -C(:,2) ...
      C(:,3)+C(:,1), -C(:,1), -C(:,3) ...
      C(:,1)+C(:,2), -C(:,2), -C(:,1)];

    assert(all(size(LI)==size(LJ)));
    assert(all(size(LI)==size(LV)));
    % Throw contribution at each edge
    L = sparse(LI,LJ,LV,3*m,n);

    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    % Map duplicate edges to first instance
    [E,~,EMAP] = unique(sort(allE,2),'rows');

    L = sparse(EMAP,F2E(:),1,size(E,1),3*m) * L;

    % Q: What's going on for boundary edges?
    % A: Interior edges are integrating around butterly. Boundary edges are only
    % integrating around half what would be their butterfly. They "should" also
    % integrate along themselves to enclose an area. Since they don't they amount
    % to computing minus the normal derivative.
  case 4
    
    warning(['This is wangyu implemetation, up to a factor with alec. Lf_alec = 2* Lf_wangyu; ']);  
      
    T = F;
    % number of vertices
    n = size(V,1);
    % number of faces
    m = size(T,1);
    % Compute cotangents
    C = cotangent(V,T);
    % Map each face-edge to a unique edge
    T2E = reshape(1:4*m,m,4);
    % Assemble entries
    LI = [ T2E(:,[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]) ];
    LJ = [   T(:,[1 2 3 4 2 3 4 1 3 4 1 2 4 1 2 3]) ];
    % allE 1st    2 2 2 
    % allE 2nd    3 3 3
    LV = [ ...
      C(:,3)+C(:,2)+C(:,4), -C(:,3), -C(:,2), -C(:,4), ...
      C(:,1)+C(:,5)+C(:,3), -C(:,1), -C(:,5), -C(:,3), ...
      C(:,6)+C(:,2)+C(:,1), -C(:,6), -C(:,2), -C(:,1), ...
      C(:,4)+C(:,5)+C(:,6), -C(:,4), -C(:,5), -C(:,6)];

    assert(all(size(LI)==size(LJ)));
    assert(all(size(LI)==size(LV)));
    % Throw contribution at each edge
    L = sparse(LI,LJ,LV,4*m,n);

    %allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    allF = [ ...
      T(:,2) T(:,4) T(:,3); ...
      T(:,1) T(:,3) T(:,4); ...
      T(:,1) T(:,4) T(:,2); ...
      T(:,1) T(:,2) T(:,3); ...
      ];
    % Map duplicate edges to first instance
    [E,~,EMAP] = unique(sort(allF,2),'rows');
    % we have allF = E(EMAP,:);
    L = sparse(EMAP,T2E(:),1,size(E,1),4*m) * L;  
      
      
%     [Lcr,E] = crouzeix_raviart_cotmatrix(V,F); 
%     A = sparse(E(:),repmat(1:size(E,1),1,3)',1,size(V,1),size(E,1))';
%     Df = diag(sparse(sum(A,2)));
%     L = 3*Lcr*(Df\A);
    % Lv == 0.5*A'*Lf;
  end

end
