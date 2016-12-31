function [B] = OneDimBiFacetLaplacian(V,Line)
    % the two implmentation are identical
    if(false)
        [M] = OneDimMass(V,Line);
        [Lz] = OneDimSquareFacetLaplacian(V,Line);
        B = Lz'*(M\Lz);
    else
        assert(size(V,2)==1);
        n = size(V,1);
        dx = V(Line(1,2),:) - V(Line(1,1),:);
        if(true)
          B = sparse(n,n);
          for i=1:n-2
              B(i:i+2,i:i+2) = B(i:i+2,i:i+2) + [1,-2,1;-2,4,-2;1,-2,1];
          end
        else
          B = sparse(n-2,n);
          for i=1:n-2
              B(i,i:i+2) = B(i,i:i+2) + [1,-2,1];
          end
          B = B'*B;
        end
        %for i=1:size(Lines,1)-1
        %  B(i:i+1,i:i+1) = B(i:i+1,i:i+1) + [1,-2,1;-2,4,-2;1,-2,1];
        %end     
        B = B/(dx*dx*dx);
    end
end