function [tv] = total_volume(V,F)

wedge = cross(V(F(:,3),:)-V(F(:,2),:),V(F(:,1),:)-V(F(:,2),:),2);

prod = 1/3 * sum( wedge .* ( V(F(:,1),:)+ V(F(:,2),:)+V(F(:,3),:) ), 2);

tv = sum(prod);