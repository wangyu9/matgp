function [Lz] = OneDimSquareFacetLaplacian(V,Lines)
    [Lez,E] = OneDimRectFacetLaplacian(V,Lines);
    n = size(V,1);
    Lz = sparse(n,n);
    Lz(E,:) = Lez;
end

