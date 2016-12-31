function [V,T,F] = my_readMESH( filename )
  % readMESH reads an MESH file with vertex/face/tet information
  %
  % [V,T,F] = readMESH( filename )
  %
  % Input:
  %  filename  path to .mesh file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  T  #T by 4 list of tet indices
  %  F  #F by 3 list of triangle indices
  %  E  #E by 2 list of tet indices

  fp = fopen(filename,'r');

  % First line is mandatory header
  MESHheader = eat_comments(fp,'#');
  check_field(MESHheader, 'MeshVersionFormatted 1');
%   if(strcmp(MESHheader,'MeshVersionFormatted 1')==0)
%     %TODO fix this later 
%     if(size('MeshVersionFormatted 1',2)==size(MESHheader,2)-1)
%     else
%     warning('First line should be "MeshVersionFormatted 1"...');
%     end
%   else
%     % force read line feed
%     fscanf(fp,'\n');
%   end
  
  Dimension3 = eat_comments(fp,'#');
  
  % second line is mandatory Dimension 3
  if(strcmp(Dimension3,'Dimension 3')==0)
    % tetgen likes to put the 3 on the next line
    % try to append next word hoping its a 3
    % force read line feed
    % wangyu do not care this: Dimension3 = [Dimension3 ' ' eat_comments(fp,'#')];
    if(strcmp(Dimension3,'Dimension 3')~=0)
      warning('Second line should be "Dimension 3"...');
    end
  end
  Vertices = eat_comments(fp,'#');
  % thrid line is mandatory Vertices
  if(strcmp(Vertices,'Vertices')~=0)
    warning('Third line should be "Vertices"...');
  end
  % read vertex count
  num_vertices = fscanf(fp,'%d\n',1);
  % read num_vertices many sets of vertex coordinates (x,y,z,ref)
  V = fscanf(fp,'%g',4*num_vertices);
  V = reshape(V,4,num_vertices)';
  V = V(:,1:3);

  Triangles = eat_comments(fp,'#');
  % forth non numeric is mandatory Triangles
  if(strcmp(Triangles,'Triangles')~=0)
    warning('Fourth (non-number) line should be "Triangles"...');
  end
  % read triangle count
  num_triangles = fscanf(fp,'%d\n',1);
  % read num_triangles many sets of face indices (a,b,c,ref)
  F = fscanf(fp,'%d',4*num_triangles);
  F = reshape(F,4,num_triangles)';
  F = F(:,1:3);

  Tetrahedra = eat_comments(fp,'#');
  % forth non numeric is mandatory Tetrahedra
  if(strcmp(Tetrahedra,'Tetrahedra')~=0)
    warning('Fifth (non-number) line should be "Tetrahedra"...');
  end
  % read tetrahedra count
  num_tetrahedra = fscanf(fp,'%d\n',1);
  % read num_tetrahedra many sets of tet indices (a,b,c,d,ref)
  T = fscanf(fp,'%d',5*num_tetrahedra);
  T = reshape(T,5,num_tetrahedra)';
  T = T(:,1:4);
  
  
  fclose(fp);
end

function check_field(MESHheader,namestring)
  if(strcmp(MESHheader,namestring)==0)
    %TODO fix this later 
    if(size(namestring,2)==size(MESHheader,2)-1)
    else
        error(['First line should be ',namestring,'...']);
    end
  else
    % force read line feed
    fscanf(fp,'\n');
  end
end


