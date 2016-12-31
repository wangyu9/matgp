function [W,pre]=BilaplacianCoordinates(V,TF,iC,varargin)
% V is all vertices
% F is all faces (2D) or tets (3D)
% iC stacks the indices of control vertices in orignal mesh.

% number of mesh vertices
n = size(V, 1);
% number of contorl handels
m = size(iC,1);


% varargin parser
pre = [];
M = speye(m);
solver = 'cholmod';

nvar = numel(varargin);
ii = 1;
while(ii<=nvar)
    switch(varargin{ii})
        case 'M'
            M = varargin{ii+1};
            ii = ii + 1;
        case 'edge_bilap'
            % should be handled at lower level in BilaplacianCoordinatesPre
            % edge_bilap = varargin{ii+1};
            ii = ii + 1;
        case 'solver'
            solver = varargin{ii+1};
            ii = ii + 1;
        case 'pre'
            pre = varargin{ii+1};
            ii = ii + 1;
    end
    ii = ii + 1;
end

known = iC;
unknown = find(~sparse(1,known,true,1,n));

if(isempty(pre))
   pre = BilaplacianCoordinatesPre(V,TF,iC,varargin{:}); 
end

% TODO: if(known~=pre.known); pre =... end


W = zeros(n,size(M,2));

fprintf('Bilaplacian Coorinate Linear Solver: ');
% tic

W(known,:) = M(:,:);

W(unknown,:) = linear_solver( pre.Qi(unknown,unknown), -pre.Qi(unknown,known)*M(:,:), solver);


% time = toc
% fprintf('Time per handle: ');
% time/length(known)

%assert(max(max(abs(V-W*V(iC,:))))<1e-7);% Double Check Linear Precision,
%should not use this with M.







% assert(max(iC)<=size(V,1));
% assert(size(iC,2)==1);
% 
% % number of mesh vertices
% n = size(V, 1);
% % number of contorl handels
% m = size(iC,1);
% 
% 
% % varargin parser
% M = speye(m);
% edge_bilap = true;
% solver = 'cholmod';
% 
% nvar = numel(varargin);
% ii = 1;
% while(ii<=nvar)
%     switch(varargin{ii})
%         case 'M'
%             M = varargin{ii+1};
%             ii = ii + 1;
%         case 'edge_bilap'
%             edge_bilap = varargin{ii+1};
%             ii = ii + 1;
%         case 'solver'
%             solver = varargin{ii+1};
%             ii = ii + 1;
%     end
%     ii = ii + 1;
% end
% 
% 
% known = iC;
% unknown = find(~sparse(1,known,true,1,n));
% 
% if(size(TF,2)==4)
%     fprintf('Solving over volume...\n');
%     
%     if(edge_bilap)
%         
%         % Compute *facet* Laplacian
%         [Lf,F] = facet_laplacian(V,TF);
%         % Mask for boundary *facets*
%         isb = is_boundary_facet(F,TF);
%         % Compute CR massmatrix for *facets*
%         [Mass,MF] = crouzeix_raviart_massmatrix(V,TF);
%         assert(all(F(:)==MF(:)));
%         % zap rows corresponding to boundary facets
%         Lfz = Lf;
%         Lfz(isb,:) = 0;
%         
%         if(true)
%             % bi-Laplacian
%             B = Lfz'*(Mass\Lfz);
%         else
%             %normalize Mass matrix
%             NM = Mass./max(Mass(:));
%             B = Lfz'*(NM\Lfz);
%         end
%         %[~,~,U] = svd(full(B));
%         %medit(V,T,boundary_faces(T),'Data',U(:,end));
%         % Linear function should be in null space:
%         
%         %assert(max(abs(B*V(:,1)))<1e-7) % wangyu turn off this assert
%         %assert(max(abs(B*V(:,1)))<1e-4*max(V(:,1)))
%         
%         is_linear_precise = max(max(abs((B*V)./V))) < 1e-4;
%         
%         if(is_linear_precise)%max(abs(B*V(:,1)))>1e-7)
%             warning(['The linear function assertion is not maintained, we should have']);
% %             max(abs(B*V(:,1)))
% %             fprintf(['<1e-7']);
% %             fprintf('Normalize Mass Matrix instead');
%         end
%     else
%         %regular laplacian
%         L = cotmatrix3(V,TF);
%         Mass = massmatrix3(V,TF,'barycentric');  
%         B = L*(Mass\L); 
%     end
% else
%     if(edge_bilap)
%             %  edge laplacian
%             [L,E] = edgeLaplacian2(V,TF);
%             % find boundary edges
%             IB = is_boundary_edge(E,TF);
%             % Construct mass matrix
%             [Mass,mE] = crouzeix_raviart_massmatrix(V,TF);
%             % Be sure same edges are being used.
%             assert(all(E(:)==mE(:)));
%             % "Kill" boundary edges
%             L(IB,:) = 0;
%     else
%             L = cotmatrix(V,TF);
%             Mass = massmatrix(V,TF,'voronoi');
%     end
%     % bilaplacian
%     B = L'*Mass*L;
% end
% 
% Qi = B;
% 
% W = zeros(n,size(M,2));
% 
% %M = sparse(b,1:m,1,n,m);
% %M = speye(m);%bc;
% 
% % invQuu = inv(Qi(unknown,unknown));
% %for i=1:m
%    % W(known,i) = Mass(:,i);
%    % W(unknown,i) = - Qi(unknown,unknown)\(Qi(unknown,known)*Mass(:,i));
%    
%    % W(unknown,i) = - invQuu*(Qi(unknown,known)*Mass(:,i));% this is very
%    % slow!
% %end
% fprintf('Bilaplacian Coorinate Linear Solver: ');
% tic
% 
% W(known,:) = M(:,:);
% 
% W(unknown,:) = linear_solver( Qi(unknown,unknown), -Qi(unknown,known)*M(:,:), solver);
% 
% 
% time = toc
% fprintf('Time per handle: ');
% time/length(known)
% 
% %assert(max(max(abs(V-W*V(iC,:))))<1e-7);% Double Check Linear Precision,
% %should not use this with M.
