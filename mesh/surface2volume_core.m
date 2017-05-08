function [V2,F2] = surface2volume_core(V1,V2,F)

% given a surface (un-closed mesh), output a volumetric mesh by shifting a
% along the direction vector d
% Then you can do [V3,T3,F3] = tetgen_with_argu(V2,F2,[],'');

% in case the input (V,F) mesh has the normals pointing to z direction
% when using the right hand rule, the normals of output mesh pointing 
% outwards.

if(size(V1,2)==2)
   V1 = [V1,zeros(size(V1,1),1)]; 
end

assert(size(F,2)==3); % input should be a 2-D manifold

[BPC11] = boundary_point_chain(F);

n = size(V1,1);

BPC21 = BPC11+n;
BPC12 = [BPC11(2:end,:);BPC11(1,:)]; %trick to make the chain "closed"
BPC22 = BPC12+n;
% connect_F = [BPC11,BPC21,BPC22;BPC11,BPC22,BPC12];
connect_F = [BPC11,BPC22,BPC21;BPC11,BPC12,BPC22]; % this is for the consideration of normal directions

V2 = [V1;V2];
F2 = [[F(:,2),F(:,1),F(:,3)];F+n;connect_F];

end

function [BPC] = boundary_point_chain(F)

    assert(size(F,2)==3); % input should be a 2-D manifold

    [I,C] = on_boundary(F);
    % this should match how the output of on_boundary is organized.
    E = [F(:,2),F(:,3); F(:,3),F(:,1); F(:,1),F(:,2)]; 
    blist = find(C(:)==1);
    BE = E(blist,:);
    
    A = adjacency_list(BE); % this function also works for Edges
    
    
    bp = size(BE,1);
    assert( bp==length(unique(BE(:))) );

    BPC =  zeros(bp,1);% boundary point chain (list)

    % initialization to select any point on the edge
    for i=1:1:bp
        if(length(A{1})>0)
           BPC(1) = i;
           break;
        end
    end

    for i=1:1:bp
        if(i<bp)
            assert(length(A{BPC(i)})==2);
            if(i==1)
                k = 1;% simple choose one direction, could be 2 as well of course.
                BPC(i+1) = A{BPC(i)}(k);
            else      
                for k=1:2
                    if(A{BPC(i)}(k)~=BPC(i-1)) % not point to it's source.
                        BPC(i+1) = A{BPC(i)}(k);
                        break;
                    end
                end
            end
        end
    end
    
    
    %BEC = zeros(size(BE));
end