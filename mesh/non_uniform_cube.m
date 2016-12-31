function [V,F] = non_uniform_cube()

N = 30;
% made an uniform cube
edge1 = [(0:1:N-1)/N*2.0-1.0, -ones(N,1)];
edge2 = [ones(N,1), (0:1:N-1)/N*2.0-1.0];
edge3 = [(N:-1:1)/N*2.0-1.0, ones(N,1)];
edge4 = [-ones(N,1), (N:-1:1)/N*2.0-1.0];

