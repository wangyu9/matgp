function [W] = NaturalNeighborCoordinates(V,H)

% number of handles
m = size(H,1);

% number of vertices to be intepolate
n = size(V,1);

W = sparse(n,m);
for i=1:m
    v = zeros(m,1);
    v(i) = 1;
    F = scatteredInterpolant(H,v,'natural');
    W(:,i) = F(V);
end

