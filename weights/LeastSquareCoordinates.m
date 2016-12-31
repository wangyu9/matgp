function [W] = LeastSquareCoordinates(V,F,b,bc,dim)

assert(dim==3); % only support 3D by now

% number of handles
m = length(b);%size(H,1);

% number of vertices to be intepolate
n = size(V,1);

known = b;
unknown = find(~sparse(1,known,true,1,n));

if(dim==2)
    error(['To be implemented!\n']);
else
    assert(dim==3);
    % in nature it is using a 2D laplacian for 3D mesh.
    L = cotmatrix(V,F);
    M = massmatrix(V,F,'voronoi');
    B = L'*(M\L);
end

Qi = B;
W = zeros(n,m);
M = speye(m);

tic
W(known,:) = M(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,:));
toc


