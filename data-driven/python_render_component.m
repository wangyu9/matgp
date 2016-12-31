function python_render_component(V,F,W)

params = [];

% get a free temporary prefix
if ~exist('prefix','var')
  prefix = tempprefix();
end

h5_file_path = [prefix '.h5'];

write_splocs(h5_file_path,W,V,F);

path_to_python = 'python C:\WorkSpace\Python\splocs\view_splocs.py'
    
command = [path_to_python ' ' params ' ' h5_file_path];

fprintf('%s\n',command);

[status, result] = system( command );

status
result

delete(h5_file_path);