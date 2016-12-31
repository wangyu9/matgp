function [W]=LocalizedBilaplacianCoordinates(V,TF,C,b,bc,edge_bilap,zero_region)

if(~exist('edge_bilap'))
    edge_bilap = true;
end

% V is all vertices
% F is all faces (2D) or tets (3D)
% C is all vertices on cage
% F_cage is all faces on cage
% We have C = M*C_input

% assume that C are embedde in V 

% number of mesh vertices
n = size(V, 1);

% number of contorl handels
m = size(C,1);

% % compute distance from every vertex in the mesh to every control vertex
% D = permute(sum((repmat(V,[1,1,c]) - ...
%     permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
% % use distances to determine closest mesh vertex to each control vertex
% % Cv(i) is closest vertex in V to ith control vertex in C
% 
% [minD,Cv] = min(D);
[row_b, col_b] = find(bc==1);
bb= b(row_b);% bb is list of handles in INPUT order.

[row_inv_b, col_inv_b] = find(inv(bc)==1);
% so we have b=bb(row_b);

if(size(TF,2)==4)
    fprintf('Solving over volume...\n');
    
    if(edge_bilap)
        
        % Compute *facet* Laplacian
        [Lf,F] = facet_laplacian(V,TF);
        % Mask for boundary *facets*
        isb = is_boundary_facet(F,TF);
        % Compute CR massmatrix for *facets*
        [M,MF] = crouzeix_raviart_massmatrix(V,TF);
        assert(all(F(:)==MF(:)));
        % zap rows corresponding to boundary facets
        Lfz = Lf;
        Lfz(isb,:) = 0;
        
        if(true)
            % bi-Laplacian
            B = Lfz'*(M\Lfz);
        else
            %normalize Mass matrix
            NM = M./max(M(:));
            B = Lfz'*(NM\Lfz);
        end
        %[~,~,U] = svd(full(B));
        %medit(V,T,boundary_faces(T),'Data',U(:,end));
        % Linear function should be in null space:
        % assert(max(abs(B*V(:,1)))<1e-7) % wangyu turn off this assert
        if(max(abs(B*V(:,1)))>1e-7)
            warning(['The linear function assertion is not maintained, we should have']);
            max(abs(B*V(:,1)))
            fprintf(['<1e-7']);
            fprintf('Normalize Mass Matrix instead');
        end
    else
        %regular laplacian
        L = cotmatrix3(V,TF);
        M = massmatrix3(V,TF,'barycentric');  
        B = L*(M\L); 
    end
else
    if(edge_bilap)
            %  edge laplacian
            [L,E] = edgeLaplacian2(V,TF);
            % find boundary edges
            IB = is_boundary_edge(E,TF);
            % Construct mass matrix
            [M,mE] = crouzeix_raviart_massmatrix(V,TF);
            % Be sure same edges are being used.
            assert(all(E(:)==mE(:)));
            % "Kill" boundary edges
            L(IB,:) = 0;
    else
            L = cotmatrix(V,TF);
            M = massmatrix(V,TF,'voronoi');
    end
    % bilaplacian
    B = L'*M*L;
end

% % harmonic system matrix
% Qi = L;%L*(M\L);
% Q = sparse(m*n,m*n);
% % Q is sparse matrix with Qi along diagonal
% for ii = 1:m
%     d = (ii - 1)*n + 1;
%     Q(d:(d + n-1), d:(d + n-1)) = Qi;
% end
% % linear constraints: partition of unity constraints and boundary
% % conditions
% PA = repmat(speye(n,n),1,m);
% Pb = ones(n,1);
% % boundary conditions
% BCAi = speye(n,n);
% BCAi = BCAi(b,:);
% BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
% % BCA is sparse matrix with BCAi along diagonal.
% for ii = 1:m
%     di = (ii - 1)*size(BCAi,1) + 1;
%     dj = (ii - 1)*size(BCAi,2) + 1;
%     BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
% end
% BCb = bc(:);
Qi = B;
W = zeros(n,m);
%M = sparse(b,1:m,1,n,m);
M = speye(m);%bc;
fprintf('Localized Bilaplacian Coorinate Linear Solver: ');
tic
% invQuu = inv(Qi(unknown,unknown));

for i=1:m
   known = b';
   add_known = setdiff(find(zero_region(:,row_inv_b(i))==1),known);
   %known = [known;add_known];
   unknown = find(~sparse(1,known,true,1,n)); 
    
   W(known,i) = zeros(size(known,1),1);%M(:,i);
   W(known(i),i) = 1;
   W(unknown,i) = - Qi(unknown,unknown)\(Qi(unknown,known)*W(known,i));
   
   % W(unknown,i) = - invQuu*(Qi(unknown,known)*M(:,i));% this is very
   % slow!
end
% W(known,:) = M(:,:);
% W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,:));
W = W(:,row_b);
time = toc
fprintf('Time per handle: ');
time/size(C,1)
