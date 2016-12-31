function W = ConstraintBilaplacianCoordinates(V,F,b,bc)
  % BIHARMONIC_BOUNDED Compute biharmonic bounded coordinates, using quadratic
  % optimizer
  %
  % W = biharmonic_bounded(V,F,b,bc,type,pou)
  % W = biharmonic_bounded(V,F,b,bc,type,pou,low,up)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %  R  region in which the weighs should be zero.
  %  Optional:
  %    type  type of optimizer to use {best available}:
  %      'quad'
  %      'least-squares'
  %      'conic'
  %    pou  true or false, enforce partition of unity explicitly {false}
  %    low  lower bound {0}
  %    up  upper bound {1}
  %  
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions
  %

energy_type = 'biharmonic';

%   if isempty(energy_type)
%     energy_type = 'biharmoic';
%   end

% number of vertices
n = size(V,1);
% number of handles
m = size(b,1);

%assert(m*n<=80000*100); % or it might run out of memory

switch(energy_type)
  case 'biharmonic'
      fprintf(['Computing biharmoic\n']);
      % Build discrete laplacian and mass matrices used by all handles' solves
      if(size(F,2)==4)
        fprintf('Solving over volume...\n');
        L = cotmatrix3(V,F);
        M = massmatrix3(V,F,'barycentric');        
      else
        L = cotmatrix(V,F);
        M = massmatrix(V,F,'voronoi');
      end
      if(false)
          %normalize Mass matrix
          NM = M./max(M(:))
          B = L * (NM\L);
      else
          B = L*(M\L); 
      end
  case'edgeLaplacian'
    assert(size(F,2)==3);

    fprintf(['Computing edge laplacian\n']);
    %  edge laplacian
    [L,E] = edgeLaplacian2(V,F);
    % find boundary edges
    IB = is_boundary_edge(E,F);
    % Construct mass matrix
    [M,mE] = crouzeix_raviart_massmatrix(V,F);
    % Be sure same edges are being used.
    assert(all(E(:)==mE(:)));
    % "Kill" boundary edges
    L(IB,:) = 0;
    % bilaplacian
    B = L'*M*L;
    lin = zeros(n,1);
  case 'harmonic'
      fprintf(['Computing hormonic\n']);
      % Build discrete laplacian and mass matrices used by all handles' solves
      if(size(F,2)==4)
        fprintf('Solving over volume...\n');
        L = cotmatrix3(V,F);
        %M = massmatrix3(V,F,'barycentric');        
      else
        L = cotmatrix(V,F);
        %M = massmatrix(V,F,'voronoi');
      end
      B = L;
end
  
% Enforce partition of unity as explicity constraints: solve for weights
% of all handles simultaneously

% biharmonic system matrix
Qi = B;%L*(M\L);
Q = sparse(m*n,m*n);
% Q is sparse matrix with Qi along diagonal
for ii = 1:m
    d = (ii - 1)*n + 1;
    Q(d:(d + n-1), d:(d + n-1)) = Qi;
end

% linear constraints: partition of unity constraints and boundary
% conditions
PA = repmat(speye(n,n),1,m);
Pb = ones(n,1);
% boundary conditions
BCAi = speye(n,n);
BCAi = BCAi(b,:);
BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
% BCA is sparse matrix with BCAi along diagonal.
for ii = 1:m
    di = (ii - 1)*size(BCAi,1) + 1;
    dj = (ii - 1)*size(BCAi,2) + 1;
    BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
end
BCb = bc(:);

      % added by wangyu
%       P = V(b,:);% get all handle vertices
%       wA = [];
%       for i=1:1:n
%         wA = sparse(blkdiag(wA,P'));
%       end
%       wB = V';
%       wB = wB(:);

P = V(b,:);% get all handle vertices
      
wAx = [];
for i=1:1:m
  wAx = [wAx,P(i,1)*speye(n,n)];
  %wAx = [wAx,1*speye(n,n)];
end
wAx(b,:) = [];%removing reproducing constrains on handles which have been included in partition of unity constraints 

wAy = [];
for i=1:1:m
  wAy = [wAy,P(i,2)*speye(n,n)];
  %wAy = [wAy,1*speye(n,n)];
end
wAy(b,:) = [];

wAz = [];
if(size(F,2)>=4)
  % 3D input
  for i=1:1:m
      wAz = [wAz,P(i,3)*speye(n,n)];               
  end
  wAz(b,:) = [];
end

wB = V;
wB(b,:) = []; 

wA = [wAx;wAy;wAz];
wB = wB(:);
% for 2D mesh input as 3D point this is useful:
wB = wB(1:size(wA,1));


M = [];
N = [];
% no sparsity constraint now
%{
group=cell(m,1);
               
for i=1:1:m
    group{i} = find(R(:,i)==1);
    meanGroup(i,:) = mean( V( group{i},:) );
end
        
for i=1:1:m       
    index_weight_zero = [];        
    index_weight_zero = group{i};

    M1 = speye(n,n);
    M1 = M1(index_weight_zero',:);
    N1 = zeros( size(M1,1), 1);
    % end
    M = blkdiag(M,M1);
    N = [N;N1];
end
%}      
% pair M,N is sparsity constraint
% pair wA,wB is coordinate reproducing constraint
      
tic;

PA(b,:) = [];
Pb(b,:) = [];

Aeq = [PA;BCA;wA;];%M  no need to add M N if using min_quad_with_fixed_zero
Beq = [Pb;BCb;wB;];%N
[rows, cols] = find(M==1);


known_zero = cols;% known zero column
tic
if(false)
    %Y = sparse(size(known,1),1);
    %[W,F,Lambda,Lambda_known] = min_quad_with_fixed(Q,zeros(n*m,1),known,Y,Aeq,Beq);

    %[known_part, ~] = find(BCA==1);
    %[W] = min_quad_with_fixed(Q,zeros(n*m,1),known_part,[BCb],[PA;wA;],[Pb;wB;]);

    [W] = min_quad_with_fixed_zero(Q,zeros(n*m,1),known_zero,Aeq,Beq);
    %[W] = min_quad_with_null_space(Q,zeros(n*m,1),known_zero,Aeq,Beq);
else
    % assert(size(known_zero)==0); %not support
    fprintf('SuiteSparse Solver ');
    [X] = linear_solver([Q,Aeq';Aeq,zeros(size(Aeq,1))],[zeros(size(Q,1),size(Beq,2));Beq],'default');
    W = X(1:n*m,1);
end
fprintf('Computing time:');toc
fprintf('\n');
            

% double check that all constraints are satisfied
remaining_value = Beq - Aeq*W;
fprintf('Double check for min_quad_with_fixed: remaining Values:\n');
remaining_indices = find(remaining_value>0.000001);
[remaining_indices,remaining_value(remaining_indices)]
toc

W = reshape(W,n,m); 
