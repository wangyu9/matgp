function [T_cluster] =  elements_clustering(TF,k,W,type)
% TF is tets or triangles
% W is what used to do clustering

%assert(max(max(TF))==size(V,1));
%assert(size(W,1)==size(V,1));

%if(exists('type')~=1)
%    type = 'farthest_sampling';
%end

assert(min(min(TF))==1);
if(min(min(TF))==0)
    TF = TF+1;warning('Add 1 to TF to make it 1-indexed');
end

dim = size(TF,2)-1;

W_for_TF = per_element_attribs_from_vertex(W,TF,'uniform');

switch(type)
    case 'farthest_sampling'
        B = farthest_element_sampling(W,TF,[],k,1);
        Dist = pdist2(W_for_TF,W_for_TF(B,:),'euclidean');
        [~,T_cluster] = min(Dist,[],2);
    case 'farthest_sampling_geo'
        assert(size(W,2)==3); % W should be V in this case
        B = farthest_element_sampling(W,TF,[],k,1,'euclidean');
        Dist = tet_geo_dist(W,TF,B);
        [~,T_cluster] = min(Dist,[],2);
    case 'kmeans'
        [T_cluster]=kmeans(W_for_TF,k);
    otherwise 
        error(['Unknown clustering method!\n']);
end

assert(min(T_cluster)==1);
assert(size(T_cluster,1)==size(TF,1));

end

function [D] = tet_geo_dist(V,TF,b) % same as farthest_element_sampling's
    T = [];
    F = [];
    n = size(V,1);
    d = size(TF,2);
    tf = size(TF,1);
    if(d==3)
        F = TF;
    else
        assert(d==4);
        T = TF;
    end
    Dv = zeros(n,size(b,1)); % dv is the vertex-tet distance
    D = zeros(tf,size(b,1)); % D is the tet-tet distance
    for i=1:1:d
        [Dv] = Dv + geo_dist(V,T,F,TF(b,i),'euclidean');
    end
    Dv = Dv / d;
    
    for i=1:1:d
       D = D +  Dv(TF(:,i),:);
    end
    D = D / d;
end