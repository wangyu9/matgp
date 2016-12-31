function [D_naive,dist_pre] = heat_geodesic_naive(V,F,b,varargin)

coeff_time = 10;

%%
nvar = numel(varargin);

ii = 1;
while(ii<=nvar)
   switch(varargin{ii})
       case 'time'
           t = varargin{ii+1};
           ii = ii + 1;
       case 'coeff_time'
           coeff_time = varargin{ii+1};
           ii = ii + 1;
       otherwise
           error('Unknown Parameter!\n');
   end
   ii = ii + 1;    
end

%%
h = mean(mean(edge_lengths(V,F))); % average edge length, warning that the internal edges are counted twice.
t = coeff_time * h^2;
%%
n = size(V,1);
m = length(b);
L = cotmatrix(V,F);
M = massmatrix(V,F,'barycentric');
rhs = sparse(b,1:m,1,n,m);
%%
%A = M\L*t;
A = L*t;
%expLtb = exp(A); % no minus due to L is negative definite already.
%expLtb = expLtb(:,b);
expLtb = (M-t*L)\M(:,b);
% first order approximation
%expLtb = (speye(n)-A)\rhs;
%c1 = 1; c2 = -1/2; c3 = 1/6;
%expLtb = (speye(n)-c1*A-c2*A^2)\rhs;
%expLtb = (speye(n)-c1*A-c2*A^2-c3*A^3)\rhs;

%D_naive = full(sqrt(-4*t*log( expLtb )));
D_naive = full(-0.5*sqrt(t)*log( expLtb )); % important: should use (5) in Keenan's paper here.

%%
dist_pre = [];% TODO
%%
D_naive = D_naive - min(D_naive);