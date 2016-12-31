function [] = point_based_skinning(V,T,F,W,H)


params = [];

% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
end

command_file_path = [prefix '.command_line.txt'];
mesh_file_path = [prefix '.obj'];
handle_file_path = [prefix 'handle.dmat'];
%weight_file_path = [prefix 'weight.dmat'];
weight_file_path = [prefix 'weight.h5'];

fh = fopen(command_file_path,'w');
%fprintf(fh,'Viewer\t\tload_mesh_from_file\t\t\t%s',mesh_file_path);
%fprintf(fh,'Viewer\t\tload_points_set_from_file\t\t\t%s',handle_file_path);
%fprintf(fh,'Viewer\t\tload_weights_from_file\t\t\t%s',weight_file_path);

add_command_line(fh,'Viewer','load_mesh_from_file',mesh_file_path);
add_command_line(fh,'Handle','load_points_set_from_file',handle_file_path);
add_command_line(fh,'Handle','add_some_purepoint_handle_from_PS');
add_command_line(fh,'Handle','clear_points_set');
add_command_line(fh,'DeformSkinning','load_weights_from_file',weight_file_path);
add_command_line(fh,'DeformPhys','set_FPS_rotation_cluster','20');
add_command_line(fh,'Handle','set_active_handles');
add_command_line(fh,'DeformPhys','start_solver');
% add_command_line(fh,);
% add_command_line(fh,);

% example command file
% Viewer			load_mesh_from_file				woody.obj
% Handle			load_points_set_from_file		H9.dmat
% Handle			add_some_purepoint_handle_from_PS
% Handle			clear_points_set
% DeformSkinning	load_weights_from_file			W9.dmat
% DeformPhys		set_FPS_rotation_cluster		20
% Handle			set_active_handles
% DeformPhys		start_solver

fclose(fh);

writeOBJ(mesh_file_path,V,F);
%writeMESH(mesh_file_path,V,T,F);

%writeDMAT(weight_file_path,W);
writeH5MAT(weight_file_path,W);

writeDMAT(handle_file_path,H);


path_to_viewer = '"C:\WorkSpace\Visual Studio 2013\PBS\viewer\msvc\x64\Release\msvc-2015-6-30.exe"';
    
command = [path_to_viewer ' ' params ' ' command_file_path];

fprintf('%s\n',command);

[status, result] = system( command );

status
result

delete(mesh_file_path);
delete(weight_file_path);
delete(handle_file_path);
delete(command_file_path);

function [] = add_command_line(fh,field_name,fun_name,varargin)
if(numel(varargin)>=1)
   if(strcmp(varargin{1},'num'))
       para = varargin{2};
       fprintf(fh,'%s\t\t%s\t\t\t%s\n',field_name,fun_name,[para]);
   else
       para = varargin{1};
       fprintf(fh,'%s\t\t%s\t\t\t%s\n',field_name,fun_name,[para]);
   end
else
   fprintf(fh,'%s\t\t%s\n',field_name,fun_name);
end
