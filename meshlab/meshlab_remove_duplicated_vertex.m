function [V2,F2] = meshlab_remove_duplicated_vertex(V,F)

[V2, F2] = meshlab_mlx(V,F,[path_to_meshlab_mlx,'remove_duplicated_vertex']);