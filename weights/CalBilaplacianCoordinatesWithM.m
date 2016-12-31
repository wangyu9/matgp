function [W_bc,C] = CalBilaplacianCoordinatesWithM(V,T,F,iC,M,edge_bilap,solver)

assert(size(iC,2)==1);

if(isempty(T))
    TF = F;
else
    TF = T;
end

[W_bc]=BilaplacianCoordinatesWithM(V,TF,iC,M,edge_bilap,solver);
C = V(iC,:);
