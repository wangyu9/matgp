function [D] = spectral_distance(V,F,b,varargin)
%%
% Biharmonic Distance: http://www.cs.princeton.edu/~funk/biharmonic.pdf
% Commute Time Distance
% DIffusion Distance

% Also see heat_divergence_distance.m, a similar function.

%%
k = 20; % number of eigenvalues used.
type = 'diffusion'; % 'diffusion', 'harmonic'=='commute', 'biharmonic'
diff_time_coeff = 1; % only used for diffusion distance
weighted_laplacian = true;

%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'k'
            k = varargin{ii+1};
            ii = ii + 1;
        case 'diff_time_coeff'
            diff_time_coeff = varargin{ii+1};
            ii = ii + 1;
        case 'weighted_laplacian'
            weighted_laplacian = varargin{ii+1};
            ii = ii + 1;
        case 'type'
            type = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end


%%
n = size(V,1);
[EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_zero_eigen',true);

%% 
% the only difference among the three distances is that they have different
% functions f(lambda), where lambda is a eigenvalue of the
% negative Laplacian operator (-L), lambda is positive.
switch(type)
    case 'diffusion'
        % require a diffusion time parameter. 
        h = mean(mean(edge_lengths(V,F))); % average edge length, warning that the internal edges are counted twice.
        t = diff_time_coeff * h^2;
        f = @(x) exp(-2*t*x); % important: x is eigenvalue that should be positive!
    case 'harmonic'
        f = @(x) 1./x;
    case 'commute'
        f = @(x) 1./x; % same as harmonic case
    case 'biharmonic'
        f = @(x) 1./x.^2;
    otherwise
        error('Unknown type of spectral distance');
end

%%
D = zeros(n,length(b));

for ii=1:k
    % important, ED(ii,ii) is the eigenvalue of negative positive Laplacian
    % L, thus -ED(ii,ii) is positive and corresponds to what is define as
    % lambda in Yaron's Biharmonic distance paper.
    
    D = D + f( -ED(ii,ii) ) * bsxfun(@minus,EV(:,ii),EV(b,ii)') .^ 2;   
end
% do not forget this!
D = sqrt(D); 
