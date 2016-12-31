function [] = write_per_element_handle_struct(fname,V,TF,num,only_small_element)

if(only_small_element)
    [I] = farthest_small_element_sampling(V,TF,[],num,1);
else
    [I] = farthest_element_sampling(V,TF,[],num,1);
end
%%
m = size(I,1);
bb = TF(I,:)';
bb = bb(:);
H = V(bb,:);
%%
%writeDMAT('Htest',H); 
GG = cell(3,1);
GG{2} = cell(0,1);
GG{3} = cell(0,1);
GG{1} = cell(m,1);
for i=1:m
   GG{1}{i} = (1:4)' + 4*(i-1); 
end

writeHANDLE(fname,V(bb,:),[],[],[],[],[],[],GG);