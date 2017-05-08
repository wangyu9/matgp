function [ZS] = zeta_signature(VV,DD,Ss,varargin)

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
    
    
    nt = length(Ss);
    ZS = zeros(n,nt);
    
    for i=1:nt
       s = Ss(i);
       % this is the right implementation:
       ZS(:,i) = VV(b,2:k).^2 * abs(1./diag(DD(2:k,2:k))).^s;
    end
    
    