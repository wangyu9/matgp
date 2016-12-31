function [C,F_C,M] = mesh_upsample(C_input,F_C_input)

M = [];

% Upsample cages, upsample each triangle
C = C_input;
F_C = F_C_input;

for i=1:1:2
    %[C,F_C] = upsample(C,F_C);
    [C,F_C,M_temp] = upsample_with_faces_index(C,F_C);
    if(i==1)
        M = M_temp;
    else
        M = M_temp*M;
    end
end

assert(max(max(C-M*C_input))<0.001);

% tsurf(F_C,C)
% axis equal;