function [W] = inverse_distance_weights(V,C,P,BE)


n = size(V,1);
p = size(P,1);
be = size(BE,1);
m = p + be;

D = zeros(n,m);

D(:,1:p) = pdist2(V,C(P,:));

K = 10; %number of samples.
for ii=1:K
    a = (ii-1)/(K-1);
    D(:,p+1:p+be) = D(:,p+1:p+be) + pdist2(V, C(BE(:,1),:)*a+C(BE(:,2),:)*(1-a) )/K;
end


D = D + 1e-9;% to prevent divided by zero
W = 1./D;

W = bsxfun(@rdivide,W,sum(W,2));
