function [V,T,F] = tetgen_with_argu(SV,SF,IV,argu)
  % TETGEN
  % [V,T,F] = tetgen(SV,SF,IV)
  %
  % Call tetgen to construct a tetrahedral volume mesh with in a given triangle
  % mesh with optional internal contrained vertices.
  %
  % Inputs:
  %   SV  list of surface vertex positions of exterior mesh, # vertices by 3
  %   SF  list of surface face indices of exterior triangle mesh, # faces by 3
  %   IV  list of internal vertex positions, # internal vertice by 3
  %   allow_resampling  allow resampling on the surface given [false]
  % Outputs:
  %   V  list of tetrahedra vertices
  %   T  list of tetrahedra indices
  %   F  list of faces of 3D volume mesh
  %

  % determine if internal constraint vertices are present
  internal_constraints = false;
  if(exist('IV','var') && prod(size(IV))>0)
    internal_constraints = true;
    assert(max(size(setdiff(SV,IV,'rows')) == size(SV)))
  end

  % get a temporary file name prefix
  prefix = tempname;
  %off_filename = [prefix '.off'];
  %writeOFF(off_filename,SV,SF);
  
  % Try to mesh with all faces included directly

  prefix = tempname;
  poly_filename = [prefix '.poly'];
  writePOLY_tetgen(poly_filename,SV,SF,[],'BoundaryMarkers',ones(size(SF,1),1));

  % if there are internal constraint vertices then print them to a .node file
  if(internal_constraints)
    inode_filename = [prefix '.a.node'];
    writeNODE(inode_filename,IV);
  end

  % graded: -q100, very-fine:-q1
  %ori: flags = '-Cp -q100 ';
%  flags = '-Cp -q50 ';
  flags = argu;

  
  if(~exist('allow_resampling','var') || ~allow_resampling)
    flags = [flags ' -Y' '-V'];
  end
  
  % flags = 'pq100Y';% wangyu
  
  
  if(internal_constraints)
    flags = [flags ' -i'];
  end
  
  % call tetgen
  path_to_tetgen = 'D:\WorkSpace2\MATLAB\bbw_demo\tetgen1.5.0\tetgen'; %wangyu
  command = [path_to_tetgen ' ' flags ' ' poly_filename];
  %fprintf(command);
  [status, result] = system(command);
  if status~=0
    error(result)
  end
  % tetgen always writes output to file:
  %   xxxx.1.ele  tetrahedra
  %   xxxx.1.node tetrahedra vertices
  %   xxxx.1.face  surface faces
  ele_filename = [prefix '.1.ele'];
  face_filename = [prefix '.1.face'];
  node_filename = [prefix '.1.node'];

  % Not sure why this isn't coming out 1-indexed
  F = readFACE(face_filename);
  F = F+1;
  % reverse faces because tetgen uses backwards order
  F = fliplr(F);
  % I guess this is 1-indexed because we're using a .off file rather than a
  % .poly file
  T = readELE(ele_filename);
  V = readNODE(node_filename);
  if min(T(:)) == 0 && max(T(:))<size(V,1)
    % make 1-indexed
    T = T + 1;
  else if min(T(:)) >= 1
    %warning('min(T) >= 1, leaving indices as is');
  end


  delete(poly_filename);
  if(internal_constraints)
    delete(inode_filename);
  end
  delete(ele_filename);
  delete(face_filename);
  delete(node_filename);

end
