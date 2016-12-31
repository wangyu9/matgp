function [M] = OneDimMass(V,Lines)
    assert(size(V,2)==1);
    n = size(V,1);
    dx = V(Lines(1,2),:) - V(Lines(1,1),:);
    M = speye(n,n) * dx;
    M(1,1) = M(1,1)*0.5;
    M(n,n) = M(n,n)*0.5;
end

