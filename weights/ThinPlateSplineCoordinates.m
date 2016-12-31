function [W,S] = ThinPlateSplineCoordinates(V,b,dim,varargin)
% see implementation details in the paper Principle Warps:
% http://user.engineering.uiowa.edu/~aip/papers/bookstein-89.pdf

kernel = 'biharmonic';
C = V(b,:);% get all handle vertices

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
   switch(varargin{ii})
       case 'kernel'
           kernel = varargin{ii+1};
           ii = ii + 1;
       case 'BC'
           BC = varargin{ii+1};
           ii = ii + 1;
       case 'C'
           C = varargin{ii+1};
           ii = ii + 1;
   end
   ii = ii + 1; 
end

if(exist('MM')~=1)
   BC = eye(size(C,1)); 
end

assert(dim==2||dim==3); % only tested on 2d and 3d yet.

% number of handles
m = length(b);

% number of vertices to be intepolate
n = size(V,1);

W = zeros(n,m);


% for i=1:m
%     v = zeros(m,1);
%     v(i) = 1;
%     % F = scatteredInterpolant(H,v,'natural');
%     if(dim==2)
%         W(:,i) = interp2(H(:,1),H(:,2),v,V(:,1),V(:,2),'spline');
%     else
%         assert(dim==3);
%         W(:,i) = interp3(H(:,1),H(:,2),H(:,3),v,V(:,1),V(:,2),V(:,3),'spline');
%     end
% end

V=V(:,1:dim);
C=C(:,1:dim);

K = construct_RBF_matrix(C,C,dim,kernel);
P = [ones(size(C,1),1),C];
L = [K,P;P',zeros(1+dim,1+dim)];
M = construct_RBF_matrix(V,C,dim,kernel);
A = [M,ones(size(V,1),1),V];

S.K = K;
S.P = P;
S.L = L;
S.M = M;
S.A = A;

W = A/L;
W = W(:,1:size(C,1))*BC;