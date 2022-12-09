function [M] = lbs_linear_matrix(V,W,dim)

n = size(V,1);
m = size(W,2);

M = zeros(n,m*(dim+1));

if dim==2
    M(:,1:dim+1:end) = bsxfun(@times,V(:,1),W);
    M(:,2:dim+1:end) = bsxfun(@times,V(:,2),W);
    M(:,3:dim+1:end) = W;
else
    assert(dim==3);
    M(:,1:dim+1:end) = bsxfun(@times,V(:,1),W);
    M(:,2:dim+1:end) = bsxfun(@times,V(:,2),W);
    M(:,3:dim+1:end) = bsxfun(@times,V(:,3),W);
    M(:,4:dim+1:end) = W;
end