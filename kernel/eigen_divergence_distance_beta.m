function D = eigen_divergence_distance_beta(V,F,b,varargin)

% Also see spectral_distance.m, a similar function.

%%
n = size(V,1);
%%
k = 100;
coeff_time = 4;
full_kernel = true;
weighted_laplacian = false; % somehow it does not work if weighted laplacian is true.
modified_laplacian = false;
%% Variables Parser

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'k'
            k = varargin{ii+1};
            ii = ii + 1;
        case 'coeff_time'
            coeff_time = varargin{ii+1};
            ii = ii + 1;
        case 'full_kernel'
            full_kernel = varargin{ii+1};
            ii = ii + 1;
        case 'weighted_laplacian'
            weighted_laplacian = varargin{ii+1};
            ii = ii + 1;
        case 'modified_laplacian'
            modified_laplacian = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

%%
h = mean(mean(edge_lengths(V,F))); % average edge length, warning that the internal edges are counted twice.
t = coeff_time * h^2;
%%
M = massmatrix(V,F,'barycentric');

if(full_kernel)
    if(false)
        assert(size(F,2)==3);
        [L,~] = facet_laplacian_pervertex(V,F);
    else
        L = cotmatrix(V,F);
    end
    EK = pinv(L);
else
    warning('Figure out why approximate kernel is not as good as full kernel');
    zero_skipped = false; % true does not work now, figure out why % true is numerically better.
    if(modified_laplacian)
        [EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_zero_eigen',true,'standard_laplacian',false);
    else
        [EV,ED] = laplacian_spectrum(V,F,k,'weighted_eigen',weighted_laplacian,'skip_zero_eigen',zero_skipped);
    end
    %EK = heat_kernel(EV,ED,t,[(1:n)'],'zero_eigen_skipped',zero_skipped);
    EK = eigen_kernel(EV,ED,[(1:n)'],'zero_eigen_skipped',zero_skipped,'p',1);
    % normalize
    %te = 1e-60;
    %EK(find(HK<te)) = te;
    %EK = bsxfun(@rdivide,EK,sum(EK,1));
end
%%
TR = @(x) exp(-2*pi.*x);
GK = @(r,s) 1/(sqrt(2*pi)*s) * exp(-(r./s).^2./2);
Ker = @(x,s) GK(TR(x),s); 
%%
sigma = 10*h;
HK = Ker(EK,sigma);
%%
% Note: we need to correct for Matlab's insistence that 0 * -Inf = NaN.
% For entropy computation, 0 * log(0) = 0. This trick works because
% min(0,NaN) = 0, and log(P) <= 0 when P <= 1.

assert(max(0,NaN)==0);

KL = @(p,q) max(0, p.*log(p./q) );
JS = @(p,q) max(0, p.*log(p./q)+q.*log(q./p) );
%AJS = @(p,q) max(0, p.*log(2*p./(q+p))+q.*log(2*q./(p+q)) ); % this kernel
% (Jensen-Shannon divergence) does not work

SD = @(p,q) (p-q).^2;

%%
%D = bsxfun(KL,expLt(:,b1)',expLt);
p = HK(:,b);
w = diag(M);
%D_sq = 2*sqrt(2*t)*full(sum( bsxfun(KL,(p)',HK) ,2));
D_sq = 2*sqrt(2*t)*full( bsxfun(SD,(p)',HK)*w );
if(~isreal(D_sq))
    warning(['The distance has imaginary part due to numerical errors. Simply ignore.']);
    D_sq = real(D_sq);
end
D_kl = sqrt(sum(D_sq,2));
%D = D_sq;
D = sqrt(D_sq);