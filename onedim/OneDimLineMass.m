function [D] = OneDimLineMass(V,Lines)
    assert(size(V,2)==1);
    n = size(V,1);
    dx = V(Lines(1,2),:) - V(Lines(1,1),:);
    D = speye(n-1,n-1) * dx;    
end

