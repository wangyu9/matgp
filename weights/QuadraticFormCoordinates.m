function [W,pre]=QuadraticFormCoordinates(Q,iC,varargin)
% V is all vertices
% F is all faces (2D) or tets (3D)
% iC stacks the indices of control vertices in orignal mesh.

% number of mesh vertices
n = size(Q, 1);
% number of contorl handels
m = size(iC,1);


% varargin parser
pre = [];
M = speye(m);
solver = 'matlab'%'cholmod';

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'M'
            M = varargin{ii+1};
            ii = ii + 1;
        case 'solver'
            solver = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

known = iC;
unknown = find(~sparse(1,known,true,1,n));

pre.Qi = Q;

% TODO: if(known~=pre.known); pre =... end

W = zeros(n,size(M,2));

fprintf('Bilaplacian Coorinate Linear Solver: ');
% tic

W(known,:) = M(:,:);

W(unknown,:) = linear_solver( pre.Qi(unknown,unknown), -pre.Qi(unknown,known)*M(:,:), solver);
