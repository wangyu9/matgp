function [xp] = project_weight(x)

x = max(x,0);
max_x = max(x);

if(max_x==0)
    xp = x;
else
    xp = x/max_x;
end