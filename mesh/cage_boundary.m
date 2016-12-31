function [B,BC,V2,T2,F2,Vc,Fc] = cage_boundary(Vm,Fm,V_Cage,F_Cage,k)

%%

if(true)
    [Vc,Fc,BC] = upsample_cage(V_Cage,F_Cage,k);
    [V2,T2,F2] = tetgen_with_argu(Vc,Fc,Vm,'-q2');
    F2 = F2 - 1;
else
    % this implementation seems incorrect 
    [Vc,Fc,BC] = upsample_cage(V_Cage,F_Cage,k);
    [V2,T2,F2] = tetgen_with_argu([Vc;Vm],[Fc;Fm+size(Vc,1)],[],'-q2');
    F2 = F2 - 1;
end
%%
if(true)
    % making use of the facts that
    % 1) tetgen always put Vc at first of the output vertices.
    % 2) V_Cage are in the first part of Vc.
    B = [1:size(V_Cage,1)]';
    % double check
    assert( max(max(V_Cage-V2(B,:)))<0.0001 );
else
    [B] = find_subset_vertex_indices(V2,V_Cage);
end











% old should remove later
% function [B,BC,V2,T2,F2] = cage_boundary(V,V_Cage,F_Cage)
% 
% V_C = V_Cage;
% F_C = F_Cage;
% 
% for i=1:1:2
%     %[C,F_C] = upsample(C,F_C);
%     [V_C,F_C,M_temp] = upsample_with_faces_index(V_C,F_C);
%     if(i==1)
%         M = M_temp;
%     else
%         M = M_temp*M;
%     end
% end
% 
% assert(max(max(V_C-M*V_Cage))<1e-6);
% 
% [V2,T2,F2] = tetgen([V_C;V],F_C,[],false);
% F2 = F2 - 1; % wangyu: there seems to be a bug of indexing in tetgen
% assert(min(min(F2))==1);
% 
% %%
% [CC,IA,IC] = intersect(V2,V_C,'rows');
% assert(size(IA,1)==size(V_C,1));
% % We have CC = V(IA,:) = V_input(IC,:) = 
% %   sparse(1:size(IA,1),IA,1,size(IA,1),size(V,1))*V = sparse(1:size(IC,1),IC,1)*V_input
% % so V_input =
% % sparse(1:size(IC,1),IC,1)\sparse(1:size(IA,1),IA,1,size(IA,1),size(V,1)) * V;
% [~,J] = find(sparse(1:size(IC,1),IC,1,size(IC,1),size(V2,1))\sparse(1:size(IA,1),IA,1,size(IA,1),size(V2,1))==1);
% % then we have V_input = V(J,:)
% assert( max(max(V_C-V2(J,:)))<1e-6 );
% 
% B = J;
% BC = M;