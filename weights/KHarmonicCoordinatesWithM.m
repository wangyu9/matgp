function [W]=KHarmonicCoordinatesWithM(V,TF,iC,M,k,facet_lap)
% V is all vertices
% F is all faces (2D) or tets (3D)
% iC stacks the indices of control vertices in orignal mesh.

assert(max(iC)<=size(V,1));
assert(size(iC,2)==1);

if(~exist('facet_lap'))
    facet_lap = false;
end

% number of mesh vertices
n = size(V, 1);

% number of contorl handels
m = size(iC,1);

known = iC;
unknown = find(~sparse(1,known,true,1,n));

[B]=KHarmonicMatrix(V,TF,k,facet_lap);

Qi = B;

W = zeros(n,size(M,2));

%M = sparse(b,1:m,1,n,m);
%M = speye(m);%bc;

% invQuu = inv(Qi(unknown,unknown));
%for i=1:m
   % W(known,i) = Mass(:,i);
   % W(unknown,i) = - Qi(unknown,unknown)\(Qi(unknown,known)*Mass(:,i));
   
   % W(unknown,i) = - invQuu*(Qi(unknown,known)*Mass(:,i));% this is very
   % slow!
%end
fprintf('Bilaplacian Coorinate Linear Solver: ');
tic

W(known,:) = M(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,:));

time = toc
fprintf('Time per handle: ');
time/length(known)

%assert(max(max(abs(V-W*V(iC,:))))<1e-7);% Double Check Linear Precision,
%should not use this with M.
