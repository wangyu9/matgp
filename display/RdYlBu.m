function h = RdYlBu(m)
% RdYlBu colormap
%
% From ColorBrewer2
% http://colorbrewer2.org/index.php?type=diverging&scheme=RdYlBu&n=7

if nargin < 1, m = size(get(gcf,'colormap'),1); end

% yy = [215, 48, 39; 
%       252, 141, 89; 
%       254, 224, 144; 
%       255, 255, 191; 
%       224, 243, 248; 
%       145, 191, 219; 
%       69, 117, 180]; 

yy = [165, 0, 38; 
      215, 48, 39; 
      244, 109, 67; 
      253, 174, 97; 
      254, 224, 144; 
      255, 255, 191; 
      224, 243, 248;
      171, 217, 233;
      116, 173, 209;
      69, 117, 180;
      49, 54, 149;]; 

x0 = 1:-2/(size(yy,1)-1):-1;

x1 = 1:-(2/(m-1)):-1;

h = interp1(x0, yy, x1, 'linear')/255;