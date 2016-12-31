function D = distance_from_vector_field(V,F,u)

legacy = false;
pre = [];
% precomputation for Poisson solve
pre.poisson = [];
  
L = cotmatrix(V,F);
%M = massmatrix(V,F,'barycentric');

  % Evaluate the vector field X
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
  if legacy
    [phi,pre.poisson] = min_quad_with_fixed( ...
      -L,div_X,gamma,zeros(numel(gamma),1),[],[],pre.poisson);
  else
    [phi,pre.poisson] = ...
      min_quad_with_fixed(-L,div_X,[],[],[],[],pre.poisson);
  end
  D = phi;
  % "Note that ???? is unique only up to an additive constant and should be
  % shifted such that the smallest distance value is zero."
  D = D - min(D(:));