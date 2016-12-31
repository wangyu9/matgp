% Pre-exist Parameters: V,F,T(optional),folder_path 

prepare_dependency;
C=V(B,:);
if( max(V(:,3)) < (min(V(:,3))+0.00001) )
    % 2D
    [W_bc,H,b,bc] = CalBilaplacianCoordinates(V,[],F,C);
else
    % 3D
    if(~exist('T'))
        % [V2,T,F2] = tetgen_with_argu(V,F,[],'-q2');
         [V2,T,F2] = tetgen(V,F,[]);
        F2 = F2 - 1; % wangyu: there seems to be a bug of indexing in tetgen
        assert(min(min(F2))==1);
        % writeMESH('temp.mesh',V,T,F);
    end
    if(~exist('V2'))
        V2 = V;
    end
    if(~exist('F2'))
        F2 = F;
    end
    [W_bc,H,b,bc] = CalBilaplacianCoordinates(V2,T,F,C);
    W_bc = W_bc(1:size(V,1),:);
end
W=W_bc;
success_flag = false;
if(size(W,1)>0)
    success_flag = true;
end
%%
% remaining_value = W*H-V;
% remaining_value = sparse(remaining_value.*(remaining_value>1e-6));
% fprintf('Remaining Values:\n'); % Double check whether all constrains are satisfied
% remaining_index = find(remaining_value>0.001);
% [remaining_index,remaining_value(remaining_index)];
