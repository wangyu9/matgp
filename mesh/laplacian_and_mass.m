function [L,Mass,B] = laplacian_and_mass(V,TF,facet_lap)

    if(~exist('facet_lap'))
        facet_lap = false;
    end

    if(facet_lap)
        %fprintf('Facet Laplacian...\n');
        if(size(TF,2)==4)
                %fprintf('Solving over volume...\n');

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

                if(false)
    %                 % bi-Laplacian
    %                 B = Lfz'*(Mass\Lfz);
    %             else
                    %normalize Mass matrix
                    Mass = Mass./max(Mass(:));
                    %B = Lfz'*(NM\Lfz);
                end
                
                if(true)% double checking
                    %[~,~,U] = svd(full(B));
                    %medit(V,T,boundary_faces(T),'Data',U(:,end));
                    % Linear function should be in null space:
                    B = Lfz'*(Mass\Lfz);
                    assert(max(abs(B*V(:,1)))<1e-7) % wangyu turn off this assert
                    if(max(abs(B*V(:,1)))>1e-7)
                        warning(['The linear function assertion is not maintained, we should have']);
                        max(abs(B*V(:,1)))
                        fprintf(['<1e-7']);
                        fprintf('Normalize Mass Matrix instead');
                    end
                end
                
                L = Lfz;
                
        else
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
                %B = L'*Mass\L;% wangyu I changed this and I believe this should be a typo. Original one: B = L'*Mass*L;
        end
        Mass = Mass./max(Mass(:));
        B = L'*(Mass\L);
    else
        if(size(TF,2)==4)
                %fprintf('Solving over volume...\n');
                %regular laplacian
                L = cotmatrix3(V,TF);
                Mass = massmatrix3(V,TF,'barycentric');  
        else
                L = cotmatrix(V,TF);
                Mass = massmatrix(V,TF,'voronoi');
        end
        % should we do normalization here?
        Mass = Mass./max(Mass(:));
        % because for standard laplacian matrix we have L'=L
        B = L*(Mass\L);
    end
end