function [V2,F2] = meshlab_mlx_from_file(V,F,path_to_mlx)

path_to_tmp = 'D:\WorkSpace\MATLAB\matgp-tmp\';

TEMP_IN_FILE = [path_to_tmp,'tmp-in.obj'];

TEMP_OUT_FILE = [path_to_tmp,'tmp-out.obj'];

writeOBJ(TEMP_IN_FILE,V,F);

command = [path_to_meshlab_server ' -i ' TEMP_IN_FILE ' -o ' TEMP_OUT_FILE ' -s ' path_to_mlx];

fprintf('%s\n',command);

[status, result] = system( command );

[V2,F2] = readOBJ(TEMP_OUT_FILE);