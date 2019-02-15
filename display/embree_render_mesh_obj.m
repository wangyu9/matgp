function [] = embree_render_mesh(V,F,varargin)


params = [];

prefix = 'D:\WorkSpace\temp\embree'; % specify an temp file prefix.

% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
end

mesh_file_path = [prefix '.obj'];
ecs_file_path = [prefix '.ecs'];
mtl_file_path = [prefix '.mtl'];

writeSCENE_obj(mesh_file_path,mtl_file_path,V,F);
[~,name,~] = fileparts(mesh_file_path);
writeECS(ecs_file_path,[name,'.obj']);

path_to_viewer = 'C:\WorkSpace\Tools\embree\Embree_v2.15.0_x64\bin\pathtracer.exe';

command = [path_to_viewer ' -c ' ecs_file_path ' --verbose 10'];


fprintf('%s\n',command);

[status, result] = system( command );

status
result

if 0
delete(mesh_file_path);
delete(ecs_file_path);
delete(mtl_file_path);
end