function L = cotmatrix_from_edges(V,F,varargin)

    l = [ ...
      sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
      sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
      sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
      ];
  
  % wangyu: sum(sum(abs(cotmatrix_from_edges(V,F)-cotmatrix(V,F))))
  % made from cotmatrix.m
    
  % L = cotmatrix(V,F)
  % L = cotmatrix(V,T)

  ss = size(F,2);
  switch ss
  case 3
    %% Could just replace everything with:

    % should change code below, so we don't need this transpose
    if(size(F,1) == 3)
      warning('F seems to be 3 by #F, it should be #F by 3');
    end
    F = F';

    % renaming indices of vertices of triangles for convenience
    i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
    % #F x 3 matrices of triangle edge vectors, named after opposite vertices
    %v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);
    % computing *unsigned* areas 
    dblA = doublearea_intrinsic(l);
    % cotangents and diagonal entries for element matrices
    cot12 = (l(:,1).^2+l(:,2).^2-l(:,3).^2)./dblA; cot23 = (l(:,2).^2+l(:,3).^2-l(:,1).^2)./dblA; cot31 = (l(:,3).^2+l(:,1).^2-l(:,2).^2)./dblA;
    % diag entries computed from the condition that rows of the matrix sum up to 1
    % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
    diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
    % indices of nonzero elements in the matrix for sparse() constructor
    i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
    j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
    % values corresponding to pairs form (i,j)
    v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3]/4; % wangyu somehow I have to divide 4 rather than 2 to match alec's one.  
    % for repeated indices (i,j) sparse automatically sums up elements, as we
    % want
    L = sparse(i,j,v,size(V,1),size(V,1));
  case 4
      assert(false) % have not been written
      
    if(size(F,1) == 4 && size(F,2) ~=4)
      warning('F seems to be 4 by #F, it should be #F by 4');
    end
    % number of mesh vertices
    n = size(V,1);
    % cotangents of dihedral angles
    C = cotangent(V,F);
    %% TODO: fix cotangent to have better accuracy so this isn't necessary
    %% Zero-out almost zeros to help sparsity
    %C(abs(C)<10*eps) = 0;
    % add to entries
    L = sparse(F(:,[2 3 1 4 4 4]),F(:,[3 1 2 1 2 3]),C,n,n);
    % add in other direction
    L = L + L';
    % diagonal is minus sum of offdiagonal entries
    L = L - diag(sum(L,2));
    %% divide by factor so that regular grid laplacian matches finite-difference
    %% laplacian in interior
    %L = L./(4+2/3*sqrt(3));
    %% multiply by factor so that matches legacy laplacian in sign and
    %% "off-by-factor-of-two-ness"
    %L = L*0.5;
    % flip sign to match cotmatix.m
    if(all(diag(L)>0))
      warning('Flipping sign of cotmatrix3, so that diag is negative');
      L = -L;
    end
  end
end
