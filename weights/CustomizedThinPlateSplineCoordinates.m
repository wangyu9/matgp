function [W] = CustomizedThinPlateSplineCoordinates(V,D,b,dim,kernel)
% see implementation details in the paper Principle Warps:
% http://user.engineering.uiowa.edu/~aip/papers/bookstein-89.pdf

assert(dim==2||dim==3); % only tested on 2d and 3d yet.

%kernel = 'biharmonic';

% number of handles
m = length(b);

% number of vertices to be intepolate
n = size(V,1);

assert(size(D,1)==n);

W = zeros(n,m);

V=V(:,1:dim);

C = V(b,:);% get all handle vertices
K = RadialBasisFunction( D(b,:), kernel );
%K = construct_RBF_matrix(C,C,dim,kernel);
P = [ones(size(C,1),1),C];
L = [K,P;P',zeros(1+dim,1+dim)];
M = RadialBasisFunction( D(:,:), kernel );
%M = construct_RBF_matrix(V,C,dim,kernel);
A = [M,ones(size(V,1),1),V];
%W = A*inv(L);
W = A/L;
W = W(:,1:size(C,1));
