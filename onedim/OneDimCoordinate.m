function W = OneDimCoordinate(Qi,V,iC,M)
  % KHARMONIC k-harmonic coordinates, "Harmonic Coordinates for Character
  % Articulation" by Joshi et al, and "An Intuitive Framework for Real-Time
  % Freeform Modeling" section 3 "Precomputed basis functions" by Botsch and
  % Kobbelt
  %
  % W = kharmonic(V,F,b,bc)
  % W = kharmonic(V,F,b,bc,k);
  %
  % Inputs:
  %  V  list of vertex positions
  %  TF  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  iC  list of boundary vertices
  %  M list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %  k  power of laplacian {k=1 %harmonic}
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights


  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(iC,2);

  known = iC;
  unknown = find(~sparse(1,known,true,1,n));
  
  %[Qi]=KHarmonicMatrix(V,TF,k,facet_lap);

  W = zeros(n,size(M,2));


fprintf('KHarmonicCoordinate Linear Solver: ');
tic

W(known,:) = M(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,:));

time = toc
fprintf('Time per handle: ');
time/length(known)

