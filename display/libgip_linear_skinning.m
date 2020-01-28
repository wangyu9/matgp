function [] = libgip_linear_skinning(folder_path,V,W,varargin)
%%

% default values if not import from varargin

azel = [];
Shading = 1.0;
dim = 3;

nvar = length(varargin);

ii=1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'view'))
       azel = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'Shading'))
       Shading = varargin{ii+1};
       ii = ii + 1;
   else
       error(['Unknown parameter!\n']);
   end
   ii = ii + 1;
end


%%
params = [];

% prefix = 'D:\WorkSpace\temp\embree';
% % get a free temporary prefix
% if ~exist('prefix','var')
%   prefix = tempprefix();
% end

prefix = tempprefix();

% mesh_file_path = [prefix 'model.mesh'];
% weights_file_path = [prefix '.weights.dmat'];
% 
% writeMESH(mesh_file_path,V,T,[]);
% writeDMAT(weights_file_path,C);

path_to_viewer = 'start /b D:\WorkSpace2\PBS\viewer\msvc\x64\Release\pbs.exe';

if ispc
    % https://stackoverflow.com/questions/8055371/how-do-i-run-two-commands-in-one-line-in-windows-cmd
    seperator = ' && ';
else
    assert(ismac)
    seperator = ' ; ';
end

folder_path = ['D:/WorkSpace/MATLAB/projects/bs/' folder_path];

writeDMAT([folder_path '/WHV.dmat'],lbs_linear_matrix(V,W,dim));

command = ['cd ' folder_path seperator  path_to_viewer ' ' folder_path '//command-matlab.txt'];

fprintf('%s\n',command);

% [status, result] = system( command );

[status,result] = system( command );

status
result

if 0
delete(mesh_file_path);
delete(color_file_path);
end

