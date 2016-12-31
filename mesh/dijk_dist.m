function [D] = dijk_dist(V,T,F,b,dist_type)
% by wangyu. This is an interface to calculate geodistance pairwisely.

if(~exist(dist_type))
    dist_type  = 'euclidean';
end

%addpath ./gaimc

% add faces edges
allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
% add tet edges
if(~isempty(T))
    tetE = [T(:,[1 2]);T(:,[2 3]);T(:,[1 3]);T(:,[1 4]);T(:,[2 4]);T(:,[3 4])];
    allE = [allE;tetE];
end
% Map duplicate edges to first instance
%[E,~,EMAP] = unique(sort(allE,2),'rows');


E = allE;
n = size(V,1);

%Dist = dist(V');% do not use this, it cal too much

A = sparse(n,n);
for k=1:size(E,1)
    i=E(k,1);
    j=E(k,2);
%     v1 = F(i,:);
%     v2 = F(j,:);
%     dist = sqrt();
    assert(i<=n);
    assert(j<=n);
    switch(dist_type)
        case 'euclidean'
            A(i,j) = sqrt(sum((V(i,:) - V(j,:)).^2));
            A(j,i) = sqrt(sum((V(i,:) - V(j,:)).^2));
        case 'connect'
            A(i,j) = 1;
            A(j,i) = 1;
    end
end

%A = spfun(@(x) x-min(min(A))+1,A); % remove the negative edges

As = convert_sparse(A);
%
% Now, we'll run Dijkstra's algorithm for every vertex and save the result
% On my 2GHz laptop, this takes 0.000485 seconds.
isequal(n,size(A,1));
D = zeros(n,length(b));

tic
for i=1:length(b)
    D(:,i) = dijkstra(As,b(i));
end
toc

% Visualize the dist disribution of Pi 
% Pi = 1;
% trisurf(F,V(:,1),V(:,2),D(:,Pi));

%{
% Let's try it without the conversion to see if we can notice the
% difference in speed.
% On my 2GHz laptop, this takes 0.001392 seconds.
D2 = zeros(n,n);
tic
for i=1:n
    D2(i,:) = dijkstra(A,i);
end
toc
%
% And just to check, let's make sure the output is the same.
isequal(D,D2)
%}