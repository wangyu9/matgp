function B = is_boundary_edge(E,F)
  % IS_BOUNDARY_EDGE Determine for each edge E if it is a boundary edge in F
  %
  % Inputs:
  %   E  #E by 2 list of edges
  %   F  #F by 3 list of triangles
  % Outputs:
  %   B  #E list bools. true iff unoriented edge occurs exactly once in F
  %     (non-manifold and non-existant edges will be false)
  %

  % http://www.mathworks.com/matlabcentral/newsreader/view_thread/165556
  allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  [~,~,EMAP]=unique(sort([E;allE],2),'rows');
  N = accumarray(EMAP,1);
  % Look of occurances of 2: one for original and another for boundary
  B = N(EMAP(1:size(E,1)))==2;

  %% Alternatively (slightly slower)
  %O = outline(F);
  %B = ismember(sort(E,2),sort(O,2),'rows');

end
