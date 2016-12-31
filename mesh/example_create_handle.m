 [C_input,F_C_input] = readOBJ('cage.obj');
 %%
 M = [];
if(true)% (tet_handle_mesh)
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
    tsurf(F_C,C)
    axis equal;
end

Cx = C(:,1);
Cy = C(:,2);
Cz = C(:,3);

% store control points in single #P by 2 list of points
C = [Cx,Cy,Cz];
%%
writeOBJ('cage-high.obj',C,F_C);
%%
writeDMAT('BC.dmat',full(M));