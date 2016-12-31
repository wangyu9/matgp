function [A_vertex] = per_vertex_attribs_from_element(A_element,TF)

assert(size(A_element,1)==size(TF,1));

n = max(max(TF));
m = size(A_element,2);

A_vertex = zeros(n,m);
connected = zeros(n,1);

for i=1:size(TF,2)
    A_vertex(TF(:,i),:) = A_vertex(TF(:,i),:) + A_element;
    connected(TF(:,i),1) = connected(TF(:,i),1) + 1;
end

for i=1:size(TF,2)
    connected( find(connected==0) ) = 1;
    A_vertex = bsxfun(@times,A_vertex,1./connected);
end