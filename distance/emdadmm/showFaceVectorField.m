function showFaceVectorField(X,T,vf)
% Shows one vector per face of a mesh, and colors the mesh by norm of the
% vector field.

Xm = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;

showDescriptor(X, T, sqrt(sum(vf.^2,2)));
axis equal; axis off;
cameratoolbar; 
hold on;

r = randperm(length(T));
m = round(length(T)/3);
ts = r(:);%1:m);
quiver3(Xm(ts,1),Xm(ts,2),Xm(ts,3),vf(ts,1),vf(ts,2),vf(ts,3),'k');