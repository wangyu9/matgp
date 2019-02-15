function [V2,F2] = meshlab_mlx(V,F,mlx_string)

[V2, F2] = meshlab_mlx_from_file(V,F,[path_to_meshlab_mlx,mlx_string,'.mlx']);

