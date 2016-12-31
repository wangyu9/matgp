function [D] = spectral_distance(V,F,b,varargin)
%%
% Biharmonic Distance: http://www.cs.princeton.edu/~funk/biharmonic.pdf
% Commute Time Distance
% DIffusion Distance

% Also see heat_divergence_distance.m, a similar function.

%%
k = 20; % number of eigenvalues used.
weighted_laplacian = true;

%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'k'  
            k = varargin{ii+1};
            ii = ii + 1;
        case 'weighted_laplacian'
            weighted_laplacian = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end


%%
n = size(V,1);
[EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_zero_eigen',true);

%% 
f = @(x) 1./x;
%%
D = zeros(n,length(b));

for ii=1:k
    D = D - f( -ED(ii,ii) ) * EV(:,ii) * EV(b,ii);  
end
D = D;% - min(D);