function [] = embree_render_mesh(V0,F0,varargin)

%%

% default values
u0 = [];
C0 = [];
filename = [];
collapse_time = [];
fast_mode = true;
source_point = [];
angle = 0;
draw_isoline = true;
draw_plane = true;
color_map = 'RdYlBu';
color_multiplier = 1;
auto_scale = true;
hold = true;
% Map of parameter names to variable names
params_to_variables = containers.Map( { 'VertexColor',  'ScaleColor',   'FileName', 'CollapseTime', 'FastMode', 'SourcePoint',  'Angle',    'DrawIsoline',  'ColorMap', 'DrawPlane',    'ColorMultiplier',  'AutoScale', 'Hold'},...
                                       {'C0',           'u0',           'filename', 'collapse_time','fast_mode','source_point', 'angle',    'draw_isoline', 'color_map','draw_plane',   'color_multiplier', 'auto_scale', 'hold'});

%% Shared Parsing Code Segment

v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous function to use assignin to this workspace 
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end

%%
VV = V0;
FF = F0;
while ~fast_mode && size(VV,1)<15000
    [VV,FF,MM] = upsample_with_faces_index(VV,FF);
    if(~isempty(u0))
        u0 = MM*u0;
    end
    if(~isempty(C0))
        C0 = MM*C0;
    end
end
V0 = VV;
F0 = FF;
%%

if(auto_scale)%min(u0)<0||max(u0)>1)
    u0s = (u0-min(u0))/(max(u0)-min(u0));
else
    u0s = u0;
end

if(size(C0,1)==1&&size(C0,2)==3)
    C0 = repmat(C0,[size(V0,1),1]);
end

if(isempty(C0)&&~isempty(u0))
    %C0 = 0.42*uint8(255*value2color(u0s,'ColorMap','default'))+1;
    %C0 = uint8(0.3*255*value2color(u0,'ColorMap','jet'))+1;
    C0 = 0.42*uint8(color_multiplier*255*value2color(u0s,'ColorMap',color_map))+1;
    %C0 = 0.42*uint8(255*value2color(u0s,'ColorMap','heat'))+1;
end

%axisangle2matrix([0 1 0],-pi/4)
%% END of Input Parsing
%% Preprocessing
V0 = V0 * axisangle2matrix([0 1 0],angle);
if ~isempty(source_point)
    source_point = source_point * axisangle2matrix([0 1 0],angle);
end    
    dis = mean(V0);
    V0 = V0 - dis;
    s = (0.5/max(max(abs(V0))));
    V0 = V0 .* s;
if ~isempty(source_point)
    source_point = source_point - dis;
    source_point = source_point .* s;
end
%%
%%
    fl = F0';
    fl = fl(:);
    V = V0(fl,:);
    C = C0(fl,:);
    F = reshape(1:3*size(F0,1),[3,size(F0,1)])';
    %%
    n = size(V,1);
    [im,UV,XY] = texture_map(F,C,n);
    im2 = imresize(im,2,'bilinear');
    %%
    imwrite(im,'D:\WorkSpace\temp\tmp.ppm');
%     imwrite(zeros(size(im)),'D:\WorkSpace\temp\zeros.ppm');
%     imwrite(ones(size(im)),'D:\WorkSpace\temp\ones.ppm');
%     imwrite(0.5*ones(size(im)),'D:\WorkSpace\temp\half.ppm');

%%
    isoline = struct();
    isoline.V = V0;
    isoline.F = F0;
    isoline.u = u0;
    isoline.draw = draw_isoline;
    isoline.source_point = source_point;
    argu = struct();
    argu.filename = filename;
    argu.collapse_time = collapse_time;
    argu.draw_plane = draw_plane;
    embree_render_mesh_core(V,F,UV,isoline,argu,hold);
end

function [im,UV,XY] = texture_map(F,C,n)
    XY = zeros(n,2);
    f = size(F,1);
    im = zeros(2*f,2,3,'uint8');
    XY(1:3:3*f-2,1) = 2*(1:f)-1;
    XY(1:3:3*f-2,2) = 1;
    im(1:2:2*f-1,1,:) = C(1:3:3*f-2,:);
    XY(2:3:3*f-1,1) = 2*(1:f)-1;
    XY(2:3:3*f-1,2) = 2;
    im(1:2:2*f-1,2,:) = C(2:3:3*f-1,:);
    XY(3:3:3*f-0,1) = 2*(1:f);
    XY(3:3:3*f-0,2) = 1;
    im(2:2:2*f-0,1,:) = C(3:3:3*f-0,:);
    im(2:2:2*f-0,2,:) = C(3:3:3*f-0,:);
    UV = XY;
    UV(:,1) = (XY(:,2) - 1);
    UV(:,2) = (XY(:,1)-1)/(2*f-1);
    %im(XY(:,1),XY(:,2),:);
end

function [im,UV,XY] = texture_map2(F,C,n)
    XY = zeros(n,2);
    f = size(F,1);
    k = 10;
    im = zeros(3*k*f,4,3,'uint8');
    
    XY(:,1) = (1:3*f)'*k + 0.5*(k-1);
    XY(:,2) = 1;
    
    tmp = repmat((1:3*f),[k,1]);
    im(:,1,:) = C(tmp(:),:);
    im(:,2,:) = C(tmp(:),:);
    im(:,3,:) = C(tmp(:),:);
    im(:,4,:) = C(tmp(:),:);
    
    UV = XY;
    UV(:,1) = 0.5;
    UV(:,2) = (XY(:,1)-1)/(3*k*f-1);
end

function [im,UV,XY] = texture_map3(F,C,n)
    XY = zeros(n,2);
    f = size(F,1);
    k = 2;
    im = zeros(k*f,k,3,'uint8');
    XY(1:3:3*f-2,1) = 2*(1:f)-1;
    XY(1:3:3*f-2,2) = 1;
    im(1:2:2*f-1,1,:) = C(1:3:3*f-2,:);
    XY(2:3:3*f-1,1) = 2*(1:f)-1;
    XY(2:3:3*f-1,2) = 2;
    im(1:2:2*f-1,2,:) = C(2:3:3*f-1,:);
    XY(3:3:3*f-0,1) = 2*(1:f);
    XY(3:3:3*f-0,2) = 1;
    im(2:2:2*f-0,1,:) = C(3:3:3*f-0,:);
    im(2:2:2*f-0,2,:) = C(3:3:3*f-0,:);
    UV = XY;
    UV(:,1) = (XY(:,2) - 1);
    UV(:,2) = (XY(:,1)-1)/(2*f-1);
    %im(XY(:,1),XY(:,2),:);
end

function [t] = texture()
t = struct();
t.map_d = [];%'ones.ppm';
t.map_Kd = 'tmp.ppm';
t.map_Ks = [];%'red.ppm';
t.map_Ns = [];%'ones.ppm';
t.map_Bump = [];%'ones.ppm';
end

function [] = embree_render_mesh_core(V,F,UV,isoline,argu,hold)

params = [];

prefix = 'D:\WorkSpace\temp\embree'; % specify an temp file prefix.

% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
end

ecs_file_path = [prefix '.ecs'];
xml_file_path = [prefix '.xml'];

writeSCENE(xml_file_path,V,F,UV,texture(),isoline,argu);
[pathstr,name,~] = fileparts(xml_file_path);

argu2 = struct();
argu2.image_size = [800,800];
writeECS(ecs_file_path,[name,'.xml'],argu2);

% path_to_viewer = 'C:\WorkSpace\Tools\embree\Embree_v2.15.0_x64\bin\pathtracer.exe';
% command = [path_to_viewer ' -c ' ecs_file_path ' --verbose 100'];

path_to_viewer = 'D:\WorkSpace\renderer\embree\build\Release\pathtracer';
if ~hold
    path_to_viewer = ['start /b ', path_to_viewer];
    % do 'start /b' so it does not block.
end
command = [path_to_viewer ' -c ' ecs_file_path ' --verbose 100 --threads 16 '];


if ~isempty(argu.collapse_time)
    command = [command, sprintf(' --collapse_time %f',argu.collapse_time)];
end

fprintf('%s\n',command);

[status, result] = system( command );

status
result

if ~isempty(argu.filename)
    [im,~] = tga_read_image(['./screenshot.tga']);
    % im = imread(); % this is not supported.
    imwrite(im,argu.filename);
    mkdir([argu.filename,'.dir']);
    copyfile(ecs_file_path,[argu.filename,'.dir/embree.ecs']);
    copyfile(xml_file_path,[argu.filename,'.dir/embree.xml']);
    copyfile('D:\WorkSpace\temp\tmp.ppm',[argu.filename,'.dir/tmp.ppm']);
    delete(['./screenshot.tga']);
end
    
if 0
delete(mesh_file_path);
delete(ecs_file_path);
delete(mtl_file_path);
end

end