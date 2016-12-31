function [W]=BilaplacianCoordinates(V,TF,iC,edge_bilap)
% V is all vertices
% F is all faces (2D) or tets (3D)
% iC stacks the indices of control vertices in orignal mesh.

assert(max(iC)<=size(V,1));
assert(size(iC,2)==1);

if(~exist('edge_bilap'))
    edge_bilap = true;
end

% number of mesh vertices
n = size(V, 1);

% number of contorl handels
m = size(iC,1);

known = iC;
unknown = find(~sparse(1,known,true,1,n));

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

Qi = B;
W = zeros(n,m);
%M = sparse(b,1:m,1,n,m);
M = speye(m);%bc;

% invQuu = inv(Qi(unknown,unknown));
%for i=1:m
   % W(known,i) = M(:,i);
   % W(unknown,i) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,i));
   
   % W(unknown,i) = - invQuu*(Qi(unknown,known)*M(:,i));% this is very
   % slow!
%end
fprintf('Bilaplacian Coorinate Linear Solver: ');
tic

W(known,:) = M(:,:);
W(unknown,:) = - Qi(unknown,unknown)\(Qi(unknown,known)*M(:,:));

time = toc
fprintf('Time per handle: ');
time/length(known)
