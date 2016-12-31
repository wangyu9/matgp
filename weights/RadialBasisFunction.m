function [U] = RadialBasisFunction(D, kernel)
% by wangyu
% construct the RBF matrix of V1 and V2, e.g. V1*V2'

if(~exist('kernel')==1)
    kernel = 'biharmonic';
end

switch(kernel)
    case 'biharmonic'
       if(true)%dim==2)
           U = -D.*D.*log(D);
           U(find(D==0)) = 0;
        else
           U = -D; % biharmonic kernel = |r|  
        end
    case 'triharmonic'
        assert(dim==3);
        U = D.^3; % triharmonic kernel = r^3
    case 'exponential'
        % Expenential Kernel
        U = exp(D.*D); % a normalizer is need, TODO.
    case 'r'
        U = D;
    case 'Gaussian'
        % Gaussian Kernel
        U=exp(-D.*D/(25*25));
    otherwise

        % Kernel 1/(r^2+constant)
        % U(i,j)=1/(temp+100);
        error(['Error: undefined kernel type!\n']);
end