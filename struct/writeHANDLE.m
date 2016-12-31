function writeHANDLE( filename, V,T,F,E,P,BE,CF,G)
  % WRITEMESH writes an MESH file with vertex/face/tet information
  %
  % writeMESH(filename,V,T,F)
  % writeMESH(filename,V,T,F,E)
  %
  % Input:
  %  filename  path to .mesh file
  %  V  #V by 3 list of vertices
  %  T  #T by 4|5 list of tet indices (additional column is color index)
  %  F  #F by 3|4 list of triangle indices (additional column is color index) 
  %  Optional:
  %    E  #E by 2|3 list of edge indices (additional column is color index) 
  %  P #P by 1 list of free point indices
  %  BE #BE by 2
  %  CE #CE by 2|3
  % varing cols, each row has num_of_vertices in BE 
  % See also readMESH
  %
  if ~exist('E','var')
    E = [];
  end
      
  %disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  % mandatory header info
  fprintf(fp,'MeshVersionFormatted 1\n');
  fprintf(fp,'Dimension 3\n');
  
  % vertices header
  fprintf(fp,'Vertices\n');
  % number of vertices
  fprintf(fp,'%d\n',size(V,1));
  % vertex positions
  fprintf(fp,'%0.15g %0.15g %0.15g 0\n',V');
  
  % triangles header
  fprintf(fp,'Triangles\n');
  % number of triangles
  fprintf(fp,'%d\n',size(F,1));
  if size(F,2) == 4
    % triangle indices + color reference
    fprintf(fp,'%d %d %d %d\n',F');
  else
    % triangle indices
    fprintf(fp,'%d %d %d 0\n',F');
  end
  
  % tetrahedra header
  fprintf(fp,'Tetrahedra\n');
  % number of tetrahedra 
  fprintf(fp,'%d\n',size(T,1));
  if size(T,2) == 5 
    % tetrahedra indices
    fprintf(fp,'%d %d %d %d %d\n',T');
  else
    % tetrahedra indices
    fprintf(fp,'%d %d %d %d 0\n',T');
  end
  
  % edge header
  fprintf(fp,'Edges\n');
  % number of tetrahedra 
  fprintf(fp,'%d\n',size(E,1));
  if size(E,2) == 3 
    % edge indices
    fprintf(fp,'%d %d %d\n',E');
  else
    % edge indices
    fprintf(fp,'%d %d 0\n',E');
  end
  
  % point header
  fprintf(fp,'Points\n');
  % number of points 
  assert( size(P,1)==0 || size(P,2)==1);
  fprintf(fp,'%d\n',size(P,1));
  % point indices
  fprintf(fp,'%d\n',P');
  
  % boneedge header
  fprintf(fp,'BoneEdges\n');
  % number of boneedges 
  assert(  isempty(BE) || size(BE,2)==2);
  fprintf(fp,'%d\n',size(BE,1));
  fprintf(fp,'%d %d\n',BE');

  % cagefacet header
  fprintf(fp,'CageFaces\n');
  % number of cagefacets 
  assert( isempty(CF) || size(CF,2)==3 );
  fprintf(fp,'%d\n',size(CF,1));
  fprintf(fp,'%d %d %d\n',CF');

  % group header
  fprintf(fp,'Groups\n');
  assert(iscell(G));
  assert(isempty(G)||size(G,2)==1);
  if(size(G,1)==3)
    GP = G{1};
    GBE = G{2};
    GCF = G{3};
  else
    GP = cell(0,0);
    GBE = cell(0,0);
    GCF = cell(0,0);
  end
  
  write_group('GroupPoints',GP,fp)
  write_group('GroupBoneEdges',GBE,fp)
  write_group('GroupCageFaces',GCF,fp)
  
  % end
  fprintf(fp,'End');
  fclose(fp);
end

function [] = write_group(namestring,G,fp)
  % header
  fprintf(fp,[namestring,'\n']);
  % number of grouppoints
  assert( isempty(G) || size(G,2)==1 );
  fprintf(fp,'%d\n',size(G,1));
  for i=1:1:size(G,1)
    assert( isempty(G{i})||size(G{i},2)==1); % G{i} should not be empty
    fprintf(fp,'%d ',G{i});
    fprintf(fp,'\n');
  end
end

