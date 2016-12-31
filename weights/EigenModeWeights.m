function W = EigenModeWeights(A,V,TF,B,BC)

n = size(V,1);
m = size(B,1);

known = B;
unknown = find(~sparse(1,known,true,1,n));

%[UU,SS,VV] = svds(A(unknown,unknown),size(unknown,1));

%[UU,SS,VV] = svds(A,n);
%Bases = UU(:,(end:-1:end+1-m)');

% [UU,SS,VV] = svds(A,m,0);
% Bases = UU(:,(1:m)');

[UU,SS,VV] = svd(full(A));
Bases = UU(:,(end:-1:end+1-m)');

if(false)
    % dense svd 
    % for 2D

    [UU,SS,VV] = svd(full(A));
    Bases = UU(:,(end:-1:end+1-m)');
    
else 
    %sparse svd 
    % for 3D
    [UU,SS,VV] = svds(A,m,0); % smallest eigen vaule first.
    Bases = UU(:,1:m);
end


C = Bases(B,:)\speye(m,m);

W = Bases*C;