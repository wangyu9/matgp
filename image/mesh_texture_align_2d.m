function [UV] = mesh_texture_align_2d(V,boxXY,boxIJ,width,height)

error('still bugs somewhere');

%[I,J,~] = find(alpha>0);
% [width, height, ~] = size(img);
% [UV] = mesh_texture_align_2d(V,[min(V);max(V)],[min(I),min(J);max(I),max(J)],width,height);


%[minI,maxI,minJ,maxJ] = boxIJ(:);%[box(1,1),box(1,2),box(2,1),box(2,2)];
%[minX,maxX,minY,maxY] = boxXY(:);
  
minI = boxIJ(1,1);
maxI = boxIJ(1,2);
minJ = boxIJ(2,1);
maxJ = boxIJ(2,2);

minX = boxXY(1,1);
maxX = boxXY(1,2);
minY = boxXY(2,1);
maxY = boxXY(2,2);

  uv = bsxfun(@times, [ minI, minJ; maxI, maxJ; minI, maxJ; maxI, minJ], 1./[width,height] );
  uv = [uv, ones(size(uv,1),1)];
  
  %uv(:,2) = 1 - uv(:,2);
  v = [ minX,maxY; maxX,minY; minX,minY; maxX,maxY; ];
  v = [v, ones(size(v,1),1)];
  
  Trans = (v'*v)\(v'*uv); % we have V*Trans = UV
  
  
  UV = [V,ones(size(V,1),1)] * Trans(:,1:2);


end