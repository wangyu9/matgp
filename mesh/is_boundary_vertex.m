function B = is_boundary_vertex(F)

  % IS_BOUNDARY_VERTEX Determine for each vertex if it is a boundary edge
  %
  %
  % Inputs:

  %   F  #F by element-size list of elements
  % Outputs:
  %   B  list of boundary vertices indices. 

  E = [];
  
  switch size(F,2)
  case 3
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    %% Alternatively (slightly slower)
    %O = outline(F);
    %B = ismember(sort(E,2),sort(O,2),'rows');
  case 4
    allE = [ ...
      F(:,2) F(:,4) F(:,3); ...
      F(:,1) F(:,3) F(:,4); ...
      F(:,1) F(:,4) F(:,2); ...
      F(:,1) F(:,2) F(:,3); ...
      ];
  end
  % http://www.mathworks.com/matlabcentral/newsreader/view_thread/165556
  [~,~,EMAP]=unique(sort([allE],2),'rows');
  
  N = accumarray(EMAP,1);
  % Look of occurances of 1
  BE = find(N(EMAP(1:size(allE,1)))==1);
  
  eallE = [allE];
  allBV = eallE(BE);
  B = unique(allBV(:));