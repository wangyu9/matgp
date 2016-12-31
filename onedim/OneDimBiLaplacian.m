function [B] = OneDimBiLaplacian(V,Line)
    [L] = OneDimLaplacian(V,Line);
    [M] = OneDimMass(V,Line);
    B = L'*(M\L);
end