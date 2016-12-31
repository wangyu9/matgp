function [W]=BilaplacianCoordinatesWithEqu(V,TF,Aeq,Beq,edge_bilap)
% V is all vertices
% F is all faces (2D) or tets (3D)

assert(size(Aeq,1)==size(Beq,1));

if(~exist('edge_bilap'))
    edge_bilap = true;
end

% number of mesh vertices
n = size(V, 1);

% number of control handels
m = size(Beq,1);

if(size(TF,2)==4)
    fprintf('Solving over volume...\n');
    
    if(edge_bilap)
        
        % Compute *facet* Laplacian
        [Lf,F] = facet_laplacian(V,TF);
        % Mask for boundary *facets*
        isb = is_boundary_facet(F,TF);
        % Compute CR massmatrix for *facets*
        [Mass,MF] = crouzeix_raviart_massmatrix(V,TF);
        assert(all(F(:)==MF(:)));
        % zap rows corresponding to boundary facets
        Lfz = Lf;
        Lfz(isb,:) = 0;
        
        if(true)
            % bi-Laplacian
            B = Lfz'*(Mass\Lfz);
        else
            %normalize Mass matrix
            NM = Mass./max(Mass(:));
            B = Lfz'*(NM\Lfz);
        end
        %[~,~,U] = svd(full(B));
        %medit(V,T,boundary_faces(T),'Data',U(:,end));
        % Linear function should be in null space:
        
        %assert(max(abs(B*V(:,1)))<1e-7) % wangyu turn off this assert
        %assert(max(abs(B*V(:,1)))<1e-4*max(V(:,1)))
        
        is_linear_precise = max(max(abs((B*V)./V))) < 1e-4;
        
        if(is_linear_precise)%max(abs(B*V(:,1)))>1e-7)
            warning(['The linear function assertion is not maintained, we should have']);
%             max(abs(B*V(:,1)))
%             fprintf(['<1e-7']);
%             fprintf('Normalize Mass Matrix instead');
        end
    else
        %regular laplacian
        L = cotmatrix3(V,TF);
        Mass = massmatrix3(V,TF,'barycentric');  
        B = L*(Mass\L); 
    end
else
    if(edge_bilap)
            %  edge laplacian
            [L,E] = edgeLaplacian2(V,TF);
            % find boundary edges
            IB = is_boundary_edge(E,TF);
            % Construct mass matrix
            [Mass,mE] = crouzeix_raviart_massmatrix(V,TF);
            % Be sure same edges are being used.
            assert(all(E(:)==mE(:)));
            % "Kill" boundary edges
            L(IB,:) = 0;
    else
            L = cotmatrix(V,TF);
            Mass = massmatrix(V,TF,'voronoi');
    end
    % bilaplacian
    B = L'*Mass*L;
end

Qi = B;

fprintf('Bilaplacian Coorinate Linear Solver: ');
tic

w = min_quad_with_fixed(Qi,[],[],[],Aeq,Beq);
W = reshape(w,[n,size(Beq,2)]);

time = toc
fprintf('Time per handle: ');
%time/length(known)
