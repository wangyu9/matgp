function [V2,F2,UV] = generate_cigar_mesh(fname)

disp = [0,20];
%%
R = 2*40;
PH = []; % handles need to be inserted inside the mesh
%PH = triC(:,1:2);
ns = 40;


theta = (1:ns)'/ns * pi;
leftC = R*[-sin(theta),-cos(theta)];
rightC = R*[sin(theta),cos(theta)];
rightC = bsxfun(@plus, rightC, [10*R,0]);
%%
V = [leftC;rightC;];
V = bsxfun(@plus,V,disp);
n = size(V,1);
E = [(1:n)',[(2:n)';1]];
%%
avg_sqr_edge_length = mean(sum((V(E(1:20,1),:)-V(E(1:20,2),:)).^2,2));

[V2,F2] = triangle([V;PH],E,[],'MaxArea',avg_sqr_edge_length/2.0,'Quality',30);

%%
UV =  bsxfun(@times, bsxfun(@minus,[V2(:,1),V2(:,2)],min(V2)), 1./max(max(V2)-min(V2)) );
%%
% 'cigar.obj'
%writeOBJ(fname,V2,F2,UV);