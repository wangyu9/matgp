function W = KHarmonicIterativeCoordinate(V,TF,iC,M,k)
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
  
  assert(k==3||k==2); % so far implemented
  
  
  if(k==3)
      fprintf('Triharmonic\n');
      [Lr,Mr,~] = laplacian_and_mass(V,TF,false);
      [Lz,~] = facet_laplacian_pervertex(V,TF);
      %B = Lz'*(Mr\Lr)*(Mr\Lz);
      % notice that A is symmetric:
      
      ze = zeros(n,n);
      A = ...
      [ze, ze, Lz';...
       ze, Lr, Mr;...
       Lz, Mr, ze];
      % solve the system A*[u;v;w;]=0 instead.
      nA = 3*n;
  else
      fprintf('Biharmonic\n');   
           [Lf,Mf,~] = laplacian_and_mass(V,TF,true);
           [Lr,Mr,~] = laplacian_and_mass(V,TF,false);
           
           [Lz,~] = facet_laplacian_pervertex(V,TF);
           
      if(false)
          ze = zeros(n,n);
           A = ...
           [ze, Lz';...
            Lz, Mr;]; 
      else
          ze = zeros(n,n);
           A = ...
           [ze, Lf';...
            Lf, Mf;]; 
      end
           
        
      assert(size(A,1)==size(A,2));
      nA = size(A,1);%2*n;
  end
  
  
  
  known = iC;
  unknown = find(~sparse(1,known,true,1,nA));
  
  uvw = zeros(nA,size(M,2));

  fprintf('KHarmonicIterative Linear Solver: ');
  tic

  if(true)
      % direct solver
      uvw(known,:) = M(:,:);
      uvw(unknown,:) = - A(unknown,unknown)\(A(unknown,known)*M(:,:));
      %uvw(unknown,:) = linear_solver(A(unknown,unknown), -(A(unknown,known)*M(:,:)), 'cholmod');
      W = uvw(1:n,:);
  elseif(false)
      rk = size(known,1);
      Aeq = sparse((1:rk)',known,ones(rk,1),rk,nA);
      Q = [A,Aeq';Aeq,zeros(rk,rk)];
      % LLT solver
      % Low = chol(Q,'lower');
      % so Low*Low' = A;
  else
      ucv = min_quad_with_fixed(A,zeros(nA,1),iC,M);
      W = uvw(1:n,:);
  end
  

      
  time = toc
  fprintf('Time per handle: ');
  time/length(known)

  