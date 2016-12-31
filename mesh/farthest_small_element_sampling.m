function [B] = farthest_small_element_sampling(V,TF,b,num,seed,varargin)

if(size(TF,2)==3)
    va = abs(doublearea(V,TF));
else
    va = abs(volume(V,TF));
end

[~,small_indices] = sort(va,'ascend');

TFU = TF(small_indices(1:num*5),:);

[B] = farthest_element_sampling(V,TFU,b,num,seed);