function [E,J] = boundary_edges(F)
  % BOUNDARY_EDGES Determine boundary edges of triangle stored in F
  %
  % E = boundary_edges(T)
  %
  % Input:
  %   F  triangle index list, m by 3, where m is the number of triangle
  % Output:
  %   E  list of boundary edges, n by 2, where n is the number of boundary
  %   edges
  %   J  list of indices into F, n by 1
  % wangyu modified from Alec's function boundary_faces
    
  assert(size(F,2)==3);
  
  % get all edges
  allE = [ ...
    F(:,1) F(:,2); ...
    F(:,2) F(:,3); ...
    F(:,3) F(:,1)];
  % sort rows so that faces are reorder in ascending order of indices
  sortedE = sort(allE,2);
  % determine uniqueness of edges
  [u,m,n] = unique(sortedE,'rows');
  % determine counts for each unique edges
  counts = accumarray(n(:), 1);
  % extract edges that only occurred once
  sorted_exteriorE = u(counts == 1,:);
  sorted_exteriorJ = m(counts == 1,:);
  % find in original faces so that ordering of indices is correct
  [I,J] = ismember(sortedE,sorted_exteriorE,'rows');
  J = mod(sorted_exteriorJ(J(I))-1,size(F,1))+1;
  E = allE(I,:);
end
