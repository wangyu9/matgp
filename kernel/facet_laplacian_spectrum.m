function [EV,ED] = laplacian_spectrum(V,F,k,varargin)

weighted_laplacian = true;
skip_first_eigen = true; % the first eigenvalue of Laplacian must be 0.
%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'weighted_eigen'
            weighted_laplacian = varargin{ii+1};
            ii = ii + 1;
        case 'skip_first_eigen'
            skip_first_eigen = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

%%

if(skip_first_eigen)
   k = k + 1; 
end

  % get cotangent matrix
  L = cotmatrix(V,F);

  % get k smallest magnitude eigenvalues and corresponding vectors
  if(weighted_laplacian)
      % This should be better, but yaron seemed to use barycentric
      % M = massmatrix(V,F);
      M = massmatrix(V,F,'barycentric');
      [EV,ED] = eigs(L,M,k,'sm');
  else
      [EV,ED] = eigs(L,k,'sm');
  end
  
  if(skip_first_eigen)
    EV = EV(:, 2:end);
    ED = ED(2:end, 2:end);
  else
    EV = EV(:, 1:end);
    ED = ED(1:end, 1:end);
  end

  % This is not exactly the same, essentially it removes the mass matrix and
  % multiplies everything by a factos of -2
  % % This also works, because of the sign change in the eigenvalues matlab
  % % reverses the output order so 0.0 is the last eigenvalue
  % [EV,ED] = eigs(-2*L,M./sum(M(:)),dim+1,'sm');
  % EV = EV(:, 1:end-1);
  % ED = ED(1:end-1, 1:end-1);

  %  divide each eigenvector by corresponding eigenvalue 
  %  divide the power by 2 first because it will appear in the denominator of
  %  distance computation *outside* the squared difference see (4) and (11)
end