function [L,E] = facet_laplacian_pervertex(V,F)
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
    
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    Im = ones(m,1);
    allCR = [ [Im*2,Im*3]; [Im*3,Im*1]; [Im*1,Im*2];];
    allOpp = [F(:,[1]);F(:,[2]);F(:,[3])];% this is the point index that is opposite to the E in the same element % this should be a 1-col vector
    allF_index = [ (1:m)';(1:m)';(1:m)'; ]; % this is the index of F where the E is living in. % this should be a 1-col vector

    % Map duplicate face-edges to first instance
    [E,~,EMAP] = unique(sort(allE,2),'rows');
    B = is_boundary_facet(E,F);
    %B(:) = false;
    
    % this is non-duplicated-inside-facets, should not use this
    % inside_facets_list = find(~B);% inside-facets (not on the boundary).
    % inside_facets = E(inside_facets_list,:);

    % this is duplicated-inside-facets
    DIF_list = find(~B(EMAP,:));
    %assert( size(find(~B),1)*2==size(DIF_list,1)); % each inside-facet should appear twice.
    DIF = allE(DIF_list,:);
    DIF_opp = allOpp(DIF_list);
    DIF_F_index = allF_index(DIF_list);
    DIF_CR = allCR(DIF_list,:); % this depends on the order in function cotangent
    
    % Assemble entries
    LI = [ DIF(:,1); DIF(:,2); DIF_opp(:); ];
    LJ = LI;
    
    % using DIF_F_index and DIF_CR(:,i) as row and column indexing
    rc = size(C,1);
    LV1 = C( DIF_F_index+rc*(DIF_CR(:,1)-1) );
    LV2 = C( DIF_F_index+rc*(DIF_CR(:,2)-1) );
    
    LV = [-LV1;-LV2;LV1+LV2 ];

    assert(all(size(LI)==size(LJ)));
    assert(all(size(LI)==size(LV)));

    L = sparse(LI,LJ,LV,n,n);

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
