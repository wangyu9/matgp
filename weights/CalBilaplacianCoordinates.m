function [W_bc,C] = CalBilaplacianCoordinates(V,T,F,iC,B,BC,edge_bilap)

assert(size(iC,2)==1);

if(isempty(T))
    TF = F;
else
    TF = T;
end

[W_bc]=BilaplacianCoordinatesWithM(V,TF,iC,M,edge_bilap);
C = V(iC,:);
