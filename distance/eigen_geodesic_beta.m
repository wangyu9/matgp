function [D_keenan,dist_pre] = heat_geodesic_keenan(V,F,b,varargin)

coeff_time_m = 10;

h = mean(mean(edge_lengths(V,F))); % average edge length, warning that the internal edges are counted twice.
t = coeff_time_m * h^2;

no_left_mass = false;
no_right_mass = false;
%%
nvar = length(varargin);
ii = 1;
while(ii<=nvar)
    if(strcmp(varargin{ii},'no_mass'))
        no_left_mass = true;
        no_right_mass = true;
    elseif(strcmp(varargin{ii},'no_left_mass'))
        no_left_mass = true;
    elseif(strcmp(varargin{ii},'no_right_mass'))
        no_right_mass = true;
    end
    ii = ii + 1;
end
%%
n = size(V,1);
L = cotmatrix(V,F);
M = massmatrix(V,F,'barycentric');

%% Heat Kernel

lM = M;
rM = M;

if(no_left_mass)
   lM = speye(n);
end

if(no_right_mass)
   rM = speye(n);
end

expLtb = (lM-t*L)\rM(:,b);

%D_naive = full(-0.5*sqrt(t)*log( expLtb )); % important: should use (5) in Keenan's paper here.
%%
k = 200;
weighted_laplacian = false;
[EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_zero_eigen',true);
f = @(x) 1./x;
EP = zeros(n,length(b));
for ii=1:k
    EP = EP - f( -ED(ii,ii) ) * EV(:,ii) * EV(b,ii);  
end

%% Keenan's heat Geodesic
u = EP;

G = grad(V,F);
Div = div(V,F);

  grad_u = reshape(G*u,size(F,1),size(V,2));
  grad_u_norm = sqrt(sum(grad_u.^2,2));
  % normalize grad_u
  normalized_grad_u = bsxfun(@rdivide,grad_u,grad_u_norm);
  % correct any zero-norm gradients
  normalized_grad_u(grad_u_norm == 0,:) = 0;
  % reverse direction
  X = -normalized_grad_u;
  
  % Solve the Poisson equation 
  % divergence of X
  div_X = Div*X(:);

    [phi] = ...
      min_quad_with_fixed(-L,div_X,[],[],[],[]);
  D_keenan = phi;
  % "Note that ???? is unique only up to an additive constant and should be
  % shifted such that the smallest distance value is zero."
  D_keenan = D_keenan - min(D_keenan(:));
%%
dist_pre = [];% TODO

