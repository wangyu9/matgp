function [] = libgip_render_mesh(V,F,varargin)
%%

% default values if not import from varargin
azel = [];
%C = [150,220,150]./255.*1.1;
C = [];% C has to be empty for later code.
u = [];
c_range = [];
nvar = length(varargin);

AO = [];
cmap = 'weights-neg';
Shading = 1.0;
VertexColor = [];
LightSource = 'default';
FaceLighting = 'gouraud';%'phong';
EdgeColor = 'none';
LightMultiplier = 1;

VP = [];
VN = [];

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
   elseif(strcmp(varargin{ii},'Quiver'))
       Quiver = varargin{ii+1};
       assert(size(Quiver,2)==6);
       VP = Quiver(:,1:3);
       VN = Quiver(:,4:6);
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


%
if(size(V,2)==2)
   warning('z axis is missing! Append zero');
   V = [V,zeros([size(V,1),1])];
end
%%
params = [];

prefix = 'D:\WorkSpace\temp\embree';
% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
end

if isempty(C)
    if isempty(u)
        u = V(:,3);
        if (max(u)-min(u))==0
            u = V(:,2);
        end
    end

    u = (u-min(u))/(max(u)-min(u));

    assert(min(u)>=0);
    assert(max(u)<=1);
    %C = value2color(u);
    ccc = colormap(my_colormap(cmap));

    C = squeeze(ind2rgb(floor(u*size(ccc,1))+1,ccc));

else
    C = 
end


mesh_file_path = [prefix '.mesh.obj'];
color_file_path = [prefix '.color.dmat'];

writeOBJ(mesh_file_path,V,F);
writeDMAT(color_file_path,C);

path_to_viewer = 'start /b D:\WorkSpace\renderer\libgip\project\build\102_visualizer_bin.exe';

command = [path_to_viewer ' ' params ' ' mesh_file_path ' ' color_file_path];

fprintf('%s\n',command);

[status, result] = system( command );

status
result

if 0
delete(mesh_file_path);
delete(color_file_path);
end

