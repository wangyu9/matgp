function [D] = commute_distance(V,F,i,dim,t)


  B = commute_embedding(V,F,dim,t);
  %B = biharmonic_embedding_yaron(V,F);
  D = sqrt(sum((repmat(B(i,:),size(B,1),1)-B(:,:)).^2,2));

  %tsurf(F,[V(:,1) V(:,2) D]);
  
end
