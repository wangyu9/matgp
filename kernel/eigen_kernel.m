function [HK] = eigen_kernel(EV,ED,b,varargin)

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
zero_eigen_skipped = false;
p = 2;
%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'zero_eigen_skipped'
            zero_eigen_skipped = varargin{ii+1};
            ii = ii + 1;
        case 'p'
            p = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

%%
assert(zero_eigen_skipped==false);
%%

n = size(EV,1);
k = size(EV,2);

assert(k==size(ED,1));
assert(k==size(ED,2));

%
if(zero_eigen_skipped)
    HK = eye(n);
    HK = HK(:,b);
else
    HK = zeros(n,length(b));
end

for ii=1:k
    % whether these is a '-' infront of ED(ii,ii) should depends on the
    % laplacian operator used to get the eigenvaludes. 
    if(zero_eigen_skipped==false && ii==1)
        continue;
    end
    HK = HK + ( - 1 ./ ED(ii,ii) ).^p * EV(:,ii) * EV(b,ii)'; 
end