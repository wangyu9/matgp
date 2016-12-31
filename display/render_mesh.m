function [t,l,h] = render_mesh(V,F,varargin)

% default values if not import from varargin
azel = [0 0];
C = [150,220,150]./255.*1.1;
w = [];

nvar = length(varargin);
ii=1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'view'))
       azel = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'FaceColor'))
       C = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ScaleColor'))
       w = varargin{ii+1};
       ii = ii + 1;
   end
   ii = ii + 1;
end


%%


if(length(w)==0)
    t = tsurf(F,V,'EdgeColor','none','FaceColor',C,'FaceLighting','phong','SpecularStrength',1.);
else
    t = tsurf(F,V,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
    
    colormap(my_colormap('weights-neg'));
    caxis([-0.2,1]);
    colorbar
    set(t,'CData',w);
    drawnow;
end

%set(t,'position',[100 100 640 360]);
view(azel); % http://www.mathworks.com/help/matlab/ref/view.html

grid off;
axis off;
camzoom(1);

whitebg([1,1,1]);


lightPos = [mean(V(:,1:2)), mean(V(:,3)) + mean((max(V(:,1:2))-min(V(:,1:2))))];

  l = [ 
       light('Position',lightPos,'Style','local');%mean(V(:,1:2)), 4*(max(V(:,1:2))-min(V(:,1:2)))
      %% light('Position',[10 -10  13],'Style','local')
       ];
    camproj('persp');
    axis equal;
    h = [];
    %h = my_add_shadow(t,l,'Ground',[0,0,-1,-50]);
%}