function [color] = value2color(u)

%% http://stackoverflow.com/questions/13885847/get-colormap-color-for-a-given-value
cmap = colormap;
i = round(1+u*(size(cmap,1)-1));
%%
color = cmap(i,:);