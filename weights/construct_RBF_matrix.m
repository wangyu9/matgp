function [U] = construct_RBF_matrix(V1, V2, dim, kernel)
% by wangyu
% construct the RBF matrix of V1 and V2, e.g. V1*V2'

if(exist('kernel')~=1)
    kernel = 'biharmonic';
end

[m,mm]=size(V1);
[n,nn]=size(V2);

U = zeros(m,n);
sq_dist = zeros(m,n);
for i=1:dim
    diff = bsxfun(@minus,V1(:,i),V2(:,i)');
    sq_dist = sq_dist + diff.*diff;
end


switch(kernel)
    case 'biharmonic'
       if(true)%dim==2)
           U = -sq_dist.*log(sq_dist);
           U(find(sq_dist==0)) = 0;
        else
           U = -sqrt(sq_dist); % biharmonic kernel = |r|  
        end
    case 'triharmonic'
        assert(dim==3);
        U = sqrt(sq_dist).^3; % triharmonic kernel = r^3
    case 'constant'
        U = ones(m,n);
    case 'r'
        U = sqrt(sq_dist);
    case 'r2'
        U = sqrt(sq_dist).^2;
    case 'r3'
        U = sqrt(sq_dist).^3;
    case 'exponential'
        % Expenential Kernel
        U = exp(sq_dist); % a normalizer is need, TODO.
    case 'logr'
        U = log(sqrt(sq_dist));
        U(find(sq_dist==0)) = 0;
    case 'r2logr'
           U = sq_dist.*log(sqrt(sq_dist));
           U(find(sq_dist==0)) = 0;
    case 'rlogr'
           U = sqrt(sq_dist).*log(sqrt(sq_dist));
           U(find(sq_dist==0)) = 0;
    case 'test1'
        U = -sq_dist.*(log(sq_dist)-1);
        U(find(sq_dist==0)) = 0;
    otherwise
        % Gaussian Kernel
        % U(i,j)=exp(-temp/(2.0*8000));
        % Kernel 1/(r^2+constant)
        % U(i,j)=1/(temp+100);
        error(['Error: undefined kernel type!\n']);
end

%{ 
% old one loop fashion, just temp code for testing
for i=1:1:m
    for j=1:1:n
        temp = V1(i,:)-V2(j,:);
        temp = temp*temp';% temp = r^2  
        switch(kernel)
            case 'biharmonic'
               if(dim==2)
                    if(temp==0)
                        U(i,j)=0;
                    else
                        U(i,j)=-temp*log(temp);
                    end
                else
                    U(i,j) = -sqrt(temp); % biharmonic kernel  
                end
            case 'triharmonic'
                assert(dim==3);
                U(i,j) = sqrt(temp)^3; % triharmonic kernel
            case 'exponential'
                % Expenential Kernel
                U(i,j)=exp(temp);
            otherwise
                % Gaussian Kernel
                % U(i,j)=exp(-temp/(2.0*8000));
                % Kernel 1/(r^2+constant)
                % U(i,j)=1/(temp+100);
                error(['Error: undefined kernel type!\n']);
        end
    end
end
%}