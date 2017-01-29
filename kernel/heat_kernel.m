function [HK] = heat_kernel(EV,ED,t,b,varargin)

% EV should be n by k, stacking the eigenvectors.
% ED should be k by k, stacking the eigenvalues in the diagonal.

% equ (3) in the paper "A Concise and Provably Informative Multi-Scale
% Signature Based on Heat Diffusion"
% http://www.lix.polytechnique.fr/~maks/papers/hks.pdf

% warning, a n by n full size heat kernel could be very slow!

% important: assume the first eigenvalue is not zero, (the zero eigenvalue is already considered when constructing HK!)
% i.e. [EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_first_eigen',true);

% Question, do columns of EV have to be normalized?

%%
skip_first = false;

%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'skip_first'
            skip_first = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

%%

n = size(EV,1);
k = size(EV,2);

assert(min(min(ED))>=0);
assert(k==size(ED,1));
assert(k==size(ED,2));

%
if(skip_first)
    error('Not implemented yet!');
    HK = eye(n);
    HK = HK(:,b);
else
    HK = zeros(n,length(b));
end

for ii=1:k
    % whether these is a '-' infront of ED(ii,ii) should depends on the
    % laplacian operator used to get the eigenvaludes. If Laplacian is 
    % negative semi-definite (i.e. eigenvaludes are negative), then there
    % should be no '-'.
    HK = HK + exp( - abs(ED(ii,ii))*t ) * EV(:,ii) * EV(b,ii)'; 
end