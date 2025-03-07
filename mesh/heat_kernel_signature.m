function [HKS] = heat_kernel_signature(VV,DD,Ts,varargin)

    n = size(VV,1);
    
    % set intial values for optional
    k = size(VV,2);
    b = 1:n;

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
    
    
    nt = length(Ts);
    HKS = zeros(n,nt);   
    for i=1:nt
       t = Ts(i);
       % HKS(:,i) = sum( VV(b,1:k).^2*exp(-t*abs(DD(1:k,1:k))), 2); this is
       % wrong since e^0 = 1, off-diagonal elements are 1!!!!
       
       % this is the right implementation:
       HKS(:,i) = VV(b,1:k).^2*exp(-t*abs(diag(DD(1:k,1:k))));  
       
    end
    
    HKS = VV(:,1:k).^2 * exp(-abs(diag(DD(1:k,1:k)))*Ts);

    