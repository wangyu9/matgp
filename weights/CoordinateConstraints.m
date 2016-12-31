function [Q,PA,Pb,BCA,BCb,wA,wB]= CoordinateConstraints(B,b,V,H)

% ususally we should have H=V(b,:)

% PA / Pb is the partition of unity constraints
% BCA / BCb is the boundary condition constraints / Lagrange Constraints
% wA / wB is the Coordinate Property / Reproducing Constraints

% Q is the mn by mn objectives.

assert(size(b,2)==1);

n = size(V,1);
m = size(b,1);
bc = eye(m);


is3D = true;

Qi = B;%L*(M\L);
Q = sparse(m*n,m*n);
% Q is sparse matrix with Qi along diagonal
for ii = 1:m
    d = (ii - 1)*n + 1;
    Q(d:(d + n-1), d:(d + n-1)) = Qi;
end

% linear constraints: partition of unity constraints and boundarynconditions
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


P = H;% get all handle vertices
      
wAx = [];
for i=1:1:m
  wAx = [wAx,P(i,1)*speye(n,n)];
end
%wAx(bb,:) = [];
      
wAy = [];
for i=1:1:m
  wAy = [wAy,P(i,2)*speye(n,n)];
end
%wAy(bb,:) = [];

wAz = [];
if(is3D)
  % 3D input
  for i=1:1:m
      wAz = [wAz,P(i,3)*speye(n,n)];               
  end
  %wAz(bb,:) = [];
end
%removing reproducing constrains on handles which have been included in partition of unity constraints 

wB = V;
%wB(bb,:) = []; 

wA = [wAx;wAy;wAz];
wB = wB(:);

