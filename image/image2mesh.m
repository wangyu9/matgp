function [V,F,UV] = image2mesh( ...
  filename, laplacian_smoothness_iterations, max_points_on_boundary, varargin)
 %changed slightly by wangyu from Alec's png2mesh
% PNG2MESH
  %
  % [V,F] = png2mesh(filename, laplacian_smoothness_iterations,
  % max_points_on_boundary)
  %
  % Mesh the non transparent part of an image. Best to have the entire shape
  % surrounded by at least one transparent pixel. Also the outer boundary and
  % boundary of holes should be roughly the same size as resampling to meet the
  % max_points_on_boundary parameter is done by edge collapsing on the boundary
  % (small boundaries will end up with too many samples and therefor the
  % triangulation will be dense there).
  % 
  % Inputs:
  %   filename  path to .png file
  %   laplacian_smoothness_iterations  Number of iterations of laplacian
  %     smoothing of the positions of boundary curve points before handing
  %     curve to Triangle
  %   max_points_on_boundary  maximum number of points in the input curve
  %     before handing to Triangle
  % Outputs:
  %   F  #faces by 3, list of face indices
  %   V  #V by 2 list of polygon vertices, note that y-coordinates will be
  %     "flipped" with respect to the image if you think of the image rows as
  %     counting from 0 at the top left corner to size(im,1) at the bottom
  %     right corner
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: bwmesh, png2poly, triangle
  %

  warning('obsolete: use bwmesh');

  basename = regexprep(filename,'\.png$','');
  
  %nvar = 1;
  ii = 1;
  while(ii<=numel(varargin))
      if(strcmp(varargin{ii},'alpha_low'))
           alpha_low = varargin{ii+1};
           ii = ii+1;
      elseif(strcmp(varargin{ii},'alpha_up'))
           alpha_up = varargin{ii+1};
           ii = ii+1;
      else
          error('Unknown parameters.\n');
      end
      ii = ii + 1;
  end
  
  % read in .png file with alpha layer
  [img,map,alpha] = imread(filename);
  if(length(alpha)==0)
      fprintf('Warning: no alpha provided in the image file, set the entire image!\n');
      alpha = 255 * ones(size(img,1),size(img,2));
  end
  
  if(exist('alpha_low')==1)
      avg_img = mean(img,3);
      alpha = alpha .* uint8(avg_img>=alpha_low);
  end
  
  if(exist('alpha_up')==1)
      avg_img = mean(img,3);
      alpha = alpha .* uint8(avg_img<=alpha_up);
      %alpha(find(img<=alpha_up)) = 0;
  end
  %
  
  % use alpha channel of png image to define poly
  [V,E,H] = image2poly(... % wangyu changed
    img,alpha,...
    laplacian_smoothness_iterations, ...
    max_points_on_boundary);
  if ~isempty(H)
    warning('Holes non-empty, but I know holes sometimes come out broken');
  end

  unr = setdiff(1:size(V,1),E(:));
  if(~isempty(unr))
    warning('Unreferenced vertices in outline...\n');
  end
  % get average squared edge length as a guess at the maximum area constraint
  % for the triangulation
  avg_sqr_edge_length = mean(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
  % triangulate the polygon
  [V,F] = triangle(V,E,H,'MaxArea',avg_sqr_edge_length/2.0,'Quality',30);
  %[V,F] = triangle(V,E,[]);
  unr = setdiff(1:size(V,1),F(:));
  if(~isempty(unr))
    warning('Removing unreferenced vertices in mesh...');
    [V,~,F] = faces_first(V,[],F);
    V=V(1:max(F(:)),:);
  end

  
  [width, height, ~] = size(img);
  
  if(false)
  
      [I,J,~] = find(alpha>0);
      uv = bsxfun(@times, [ min([I,J]); max([I,J]); min(I),max(J) ; max(I),min(J) ], 1./[width,height] );
      uv = [uv, ones(size(uv,1),1)];

      %uv(:,2) = 1 - uv(:,2);
      v = [ min(V(:,1)),max(V(:,2)); max(V(:,1)),min(V(:,2)); min(V(:,1)),min(V(:,2)); max(V(:,1)),max(V(:,2)); ];
      v = [v, ones(size(v,1),1)];

      Trans = (v'*v)\(v'*uv); % we have V*Trans = UV


      UV = [V,ones(size(V,1),1)] * Trans(:,1:2);
      %UV(:,2) = 1.0 - UV(:,2);
  
  else
     UV =  bsxfun(@times, [V(:,1),V(:,2)], 1./[height,width]);
     UV(:,2) = 1.0 - UV(:,2);
  end
  
end
