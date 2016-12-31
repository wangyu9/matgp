function [B]=KHarmonicMatrix(V,TF,k,facet_lap)
% V is all vertices
% F is all faces (2D) or tets (3D)

    if(k==4)
        fprintf('Mixed Laplacian...\n');
        if(false)
            [Lf,Mf,Bf] = laplacian_and_mass(V,TF,true);
            [~,Mr,~] = laplacian_and_mass(V,TF,false);
            % This by now is just what I hack to do experiment
            % B = Lf'*Lf*Lf'*Lf;
            % B = Bf'* (Bf);
            B = Bf'*(Mr\Bf);
        else
            [Lr,Mr,~] = laplacian_and_mass(V,TF,false);
            [Lz,~] = facet_laplacian_pervertex(V,TF);
            B = Lz'*(Mr\Lr)*(Mr\Lr)*(Mr\Lz);
        end
    elseif(k==3)
        %[Lz,~] = facet_laplacian_pervertex(V,F)
        %B = Lz'*(Mass\Lz);% L is symmetric

        fprintf('Triharmonic\n');
        [Lr,Mr,~] = laplacian_and_mass(V,TF,false);
        [Lz,~] = facet_laplacian_pervertex(V,TF);
        B = Lz'*(Mr\Lr)*(Mr\Lz);
    else
        
        assert(k==2); % so far implemented
        
        if(false)
            [~,~,B] = laplacian_and_mass(V,TF,facet_lap);
        else
            fprintf('Lz\n');
            Mass = massmatrix(V,TF,'voronoi');
            Mass = Mass./max(Mass(:));
            [Lz,~] = facet_laplacian_pervertex(V,TF);
            B = Lz'*(Mass\Lz);
        end
    end
end



