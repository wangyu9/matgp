function [G] = OneDimGradient(V,Lines)
% Gradient is #element (#Lines) by #vertices
    assert(size(V,2)==1);
    n = size(V,1);
    dx = V(Lines(1,2),:) - V(Lines(1,1),:);
    G = sparse(n-1,n);
    for i=1:n-1
      G(i,i:i+1) = G(i,i:i+1) + [-1,1];
    end
    % all values live on line 
end

