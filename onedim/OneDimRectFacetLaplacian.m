function [Lez,E] = OneDimRectFacetLaplacian(V,Lines)
    assert(size(V,2)==1);
    n = size(V,1);
    dx = V(Lines(1,2),:) - V(Lines(1,1),:);
    Lez = sparse(n-2,n);
    for i=1:n-2
      Lez(i,i:i+2) = Lez(i,i:i+2) + [1,-2,1];
    end
    % all values live on vertices not on boundary (1-D facet)
    E = [2:(n-1)]';
end

