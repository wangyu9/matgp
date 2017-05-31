function writeXML(mesh_file_path, mtl_file_path, V,F)
  % WRITEOBJ writes an OBJ file with vertex/face information
  %
  % writeOBJ(filename,V,F,UV,N)
  %
  % Input:
  %  filename  path to .obj file
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #UV by 2 list of texture coordinates
  %  TF  #TF by 3 list of corner texture indices into UV
  %  N  #N by 3 list of normals
  %  NF  #NF by 3 list of corner normal indices into N
  %

writeMTL(mtl_file_path);  
  
f = fopen( mesh_file_path, 'w' );

[~,name,~] = fileparts(mtl_file_path);

fprintf( f, ['mtllib ', name,'.mtl\n']);

%fprintf( f, 'usemtl lightred\n');
fprintf( f, 'usemtl lambert4SGtexture\n');

fprintf( f, 'o model\n');

writeOBJ_core(f,V,F,V(:,[1,2])+0.5);

fprintf( f, 'o floor\n');
fprintf( f, 'usemtl white\n');

floor_y = min(V(:,2)) - 0.002;

fprintf( f, 'v %0.17g %0.17g %0.17g\n', [ 05	floor_y -05]);
fprintf( f, 'v %0.17g %0.17g %0.17g\n', [-05	floor_y -05]);
fprintf( f, 'v %0.17g %0.17g %0.17g\n', [-05	floor_y  05]);
fprintf( f, 'v %0.17g %0.17g %0.17g\n', [ 05	floor_y  05]);
fprintf( f, 'f -4 -3 -2 -1\n');

fclose(f);


end

function writeOBJ_core(f,V,F,UV,TF,N,NF)

if size(V,2) == 2
  warning('Appending 0s as z-coordinate');
  V(:,end+1:3) = 0;
else
  assert(size(V,2) == 3);
end
fprintf( f, 'v %0.17g %0.17g %0.17g\n', V');

hasN =  exist('N','var') && ~isempty(N);
hasUV = exist('UV','var') && ~isempty(UV);

if hasUV
    switch size(UV,2)
    case 2
      fprintf( f, 'vt %0.17g %0.17g\n', UV');
    case 3
      fprintf( f, 'vt %0.17g %0.17g %0.17g\n', UV');
    end
end

if hasN
    %for k=1:size(N,1)
    %    fprintf( f, 'vn %f %f %f\n', N(k,1), N(k,2), N(k,3) );
    %end
    fprintf( f, 'vn %0.17g %0.17g %0.17g\n', N');
end

if hasUV && (~exist('TF','var') || isempty(TF))
    TF = F;
end
if hasN && (~exist('NF','var') || isempty(NF))
    NF = F;
end

for k=1:size(F,1)
    if ( (~hasN) && (~hasUV) ) || (any(TF(k,:)<=0,2) && any(NF(k,:)<=0,2))
        fmt = repmat(' %d',1,size(F,2));
        fprintf( f,['f' fmt '\n'], F(k,:));
    elseif ( hasUV && (~hasN || any(NF(k,:)<=0,2)))
        fmt = repmat(' %d/%d',1,size(F,2));
        fprintf( f, ['f' fmt '\n'], [F(k,:);TF(k,:)]);
    elseif ( (hasN) && (~hasUV || any(TF(k,:)<=0,2)))
        fmt = repmat(' %d//%d',1,size(F,2));
        fprintf( f, ['f' fmt '\n'],[F(k,:);TF(k,:)]');
    elseif ( (hasN) && (hasUV) )
        assert(all(NF(k,:)>0));
        assert(all(TF(k,:)>0));
        fmt = repmat(' %d/%d/%d',1,size(F,2));
        fprintf( f, ['f' fmt '\n'],[F(k,:);TF(k,:);NF(k,:)]);
    end
end

end