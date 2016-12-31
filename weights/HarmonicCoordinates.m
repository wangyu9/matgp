function [W]=HarmonicCoordinates(V,F,B,BC)

% V is all vertices
% F is all faces (2D) or tets (3D)
% C is all vertices on cage
% F_cage is all faces on cage
% We have C = M*C_input

% assume that C are embedde in V 

assert(size(B,2)==1);

% number of mesh vertices
n = size(V, 1);

% number of boundary vertices
% c = size(C,1);

% number of contorl handels
m = size(BC,2);

% % compute distance from every vertex in the mesh to every control vertex
% D = permute(sum((repmat(V,[1,1,c]) - ...
%     permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
% % use distances to determine closest mesh vertex to each control vertex
% % Cv(i) is closest vertex in V to ith control vertex in C
% 
% [minD,Cv] = min(D);

known = B;
unknown = find(~sparse(1,known,true,1,n));

if(size(F,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,F);
    Ma = massmatrix3(V,F,'barycentric');
else
    L = cotmatrix(V,F);
    Ma = massmatrix(V,F,'voronoi');
end

% % harmonic system matrix
% Qi = L;%L*(M\L);
% Q = sparse(m*n,m*n);
% % Q is sparse matrix with Qi along diagonal
% for ii = 1:m
%     d = (ii - 1)*n + 1;
%     Q(d:(d + n-1), d:(d + n-1)) = Qi;
% end
% % linear constraints: partition of unity constraints and boundary
% % conditions
% PA = repmat(speye(n,n),1,m);
% Pb = ones(n,1);
% % boundary conditions
% BCAi = speye(n,n);
% BCAi = BCAi(b,:);
% BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
% % BCA is sparse matrix with BCAi along diagonal.
% for ii = 1:m
%     di = (ii - 1)*size(BCAi,1) + 1;
%     dj = (ii - 1)*size(BCAi,2) + 1;
%     BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
% end
% BCb = bc(:);
Qi = L;%Ma\L;
W = zeros(n,m);
% for i=1:m
%    W(known,i) = M(:,i);
%    W(unknown,i) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,i));
% end
W(known,:) = BC;%M(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*BC);
