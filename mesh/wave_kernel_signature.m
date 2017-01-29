function [WKS] = wave_kernel_signature(VV,DD,es,varargin)

    n = size(VV,1);
    
    % set intial values for optional
    k = size(VV,2);
    b = 1:n;
    wks_variance = 6;

    % parse optional inputs
    ii = 1;
    num_ri = numel(varargin);
    while(ii <= numel(varargin))
       switch varargin{ii}
           case 'k'
               k = varargin{ii+1};
               ii = ii + 1;
           case 'b'
               b = varargin{ii+1};
               ii = ii + 1;    
           otherwise
               error(['''' varargin{ii} ''' is not a valid parameter']);
       end
       ii = ii + 1;
    end
    % end of parsing
    
    log_E=log(diag(DD(1:k,1:k)))';
    sigma=(es(2)-es(1))*wks_variance;
    
    nt = length(es);
    WKS = zeros(n,nt);
    C = zeros(1,nt);
    % 
    for i=1:nt
       e = es(i);
       WKS(:,i) = sum( VV(b,1:k).^2*diag(exp((-(e-log_E).^2)./(2*sigma.^2))), 2);
       C(i) = sum(exp((-(e-log_E).^2)./(2*sigma.^2)));
    end
    
    WKS = bsxfun(@rdivide,WKS,C);
   

