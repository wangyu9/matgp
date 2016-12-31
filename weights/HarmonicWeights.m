function [W]=HarmonicWeights(V,TF,iC,varargin)


% V is all vertices
% iC is list of boundary vertices
% BC is boundary conditions on B.
% F is all faces (2D) or tets (3D)
% iC stacks the indices of control vertices in orignal mesh.
% edge_bilap should be false most of the time, use standard laplacian is enough

assert(max(iC)<=size(V,1));
assert(size(iC,2)==1);

if(~exist('edge_bilap'))
    edge_bilap = false;
end

% number of mesh vertices
n = size(V, 1);

% number of contorl handels
m = size(iC,1);

known = iC;
unknown = find(~sparse(1,known,true,1,n));


%% Default Values 
BC = eye(m);

%% Input Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'BC'
            BC = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

%%

% note L is negative positive

if(size(TF,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,TF);
    Mass = massmatrix3(V,TF,'barycentric');
%     B0 = L*(Mass\L);
%     B1 = sparse(n,n);
%     B1(iC,iC) = B0(iC,iC);
%     B2 = L;%*(Mass\L);
%     B = B0 + B2;
    B = -L; 
else
    L = cotmatrix(V,TF);
    Mass = massmatrix(V,TF,'voronoi');
%     B0 = L*(Mass\L);
%     B1 = sparse(n,n);
%     B1(iC,iC) = B0(iC,iC);
%     B2 = L;%*(Mass\L);
%     B = B0 + B2;
    B = -L;
end

Qi = B;

W = zeros(n,size(BC,2));

fprintf('Bilaplacian Coorinate Linear Solver: ');
tic

W(known,:) = BC(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*BC(:,:));

time = toc
fprintf('Time per handle: ');
time/length(known)
