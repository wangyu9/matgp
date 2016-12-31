function [W]=BilaplacianCoordinatesWithM(V,TF,iC,M,edge_bilap,solver)

W = BilaplacianCoordinates(V,TF,iC,'M',M,'edge_bilap',edge_bilap,'solver',solver);