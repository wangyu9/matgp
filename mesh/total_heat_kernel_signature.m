function [THKS] = total_heat_kernel_signature(VV,DD,Ts,varargin)

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
           otherwise
               error(['''' varargin{ii} ''' is not a valid parameter']);
       end
       ii = ii + 1;
    end
    % end of parsing
    
    
%     nt = length(Ts);
%     HKS = zeros(n,nt);
%     
%     for i=1:nt
%        t = Ts(i);
%        HKS(:,i) = sum( VV(b,1:k).^2*exp(-t*DD(1:k,1:k)), 2);        
%     end
    
    HKS = abs(VV(:,1:k) ).^2 * exp(-abs(diag(DD(1:k,1:k)))*Ts);

    THKS = sum(HKS,1);
    