function [] = libgip_render_mesh(V,F,C)


params = [];

prefix = 'D:\WorkSpace\temp\embree';
% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
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

