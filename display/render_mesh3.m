function [t,l,h] = render_mesh3(V,F,varargin)

% default values if not import from varargin
azel = [];
C = [150,220,150]./255.*1.1;
u = [];
c_range = [];
nvar = length(varargin);

AO = [];
cmap = 'weights-neg';
Shading = 1.0;
VertexColor = [];
LightSource = 'default';
FaceLighting = 'phong';
EdgeColor = 'none';
LightMultiplier = 1;

ii=1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'view'))
       azel = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'AmbientOcclusion'))
       AO = varargin{ii+1};
       ii = ii + 1;   
   elseif(strcmp(varargin{ii},'Shading'))
       Shading = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'LightSource'))
       LightSource = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'FaceColor'))
       C = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'VertexColor'))
       VertexColor = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ScaleColor'))
       u = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ColorMap'))
       cmap = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'ColorAxis'))
       c_range = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'EdgeColor'))
       EdgeColor = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'LightMultiplier'))
       LightMultiplier = varargin{ii+1};
       ii = ii + 1;
   else
       error(['Unknown parameter!\n']);
   end
   ii = ii + 1;
end


%%

t = tsurf(F,V,'EdgeColor',EdgeColor,'FaceColor','interp','FaceLighting',FaceLighting);
t.Vertices = V*axisangle2matrix([0 0 1],pi)*axisangle2matrix([1 0 0],pi/2);
colormap(my_colormap(cmap));
if(~isempty(c_range))
   caxis([c_range(1),c_range(2)]); 
end

if(isempty(AO))
AO = ambient_occlusion(V,F,V,per_vertex_normals(V,F),1000);
end

if(~isempty(u))
   % set(t,'CData',w);
    %VertexColor = value2color(u);
    assert(min(u)>=0);
    assert(max(u)<=1);
    VertexColor = squeeze(ind2rgb(floor(u*size(colormap,1))+1,colormap));
end

if(~isempty(VertexColor))
    t.FaceVertexCData = bsxfun(@times,VertexColor,1-Shading*AO);
end

t.SpecularStrength = 0.2;
t.DiffuseStrength = 0.1;
t.AmbientStrength = 0.7;

switch(LightSource)
    case 'none'
    case 'default'

        s = 0.3*LightMultiplier;

        l = [ 
            light('Position',[+1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
          ];
end
   

h = [];
    %h = my_add_shadow(t,l,'Ground',[0,0,-1,-50]);



camproj('persp');
axis equal;
grid off;
if(~isempty(azel))
view(azel); % http://www.mathworks.com/help/matlab/ref/view.html
end
axis off;
%camzoom(1);

set(gcf,'color','w');

drawnow;
