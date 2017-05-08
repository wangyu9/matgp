function [color] = value2color(u,varargin)


% default values
cm = 'default';
% Map of parameter names to variable names
params_to_variables = containers.Map( {'ColorMap'}, {'cm'});

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

%% http://stackoverflow.com/questions/13885847/get-colormap-color-for-a-given-value
colormap(cm);
cmap = colormap;
i = round(1+u*(size(cmap,1)-1));
%%
close();
%%
color = cmap(i,:);