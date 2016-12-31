function [L] = OneDimLaplacian(V,Lines)
    [D] = OneDimLineMass(V,Lines);
    [G] = OneDimGradient(V,Lines);
    L = G'*D*G;
end

