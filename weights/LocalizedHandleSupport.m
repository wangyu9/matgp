function [LHS] = LocalizedHandleSupport(V,TF,b,k)

% LHS is a matrix n by m such that is 1 if the weight is on, 0 if off. 

assert(size(b,2)==1);

if(size(TF,2)==4)
    T = TF;
    F = [];
else
    T = [];
    F = TF;
end

if(true)
    % Geodesic Distance
    [D] = geo_dist(V,T,F,b,'euclidean');% add T here since it seems that tetgen won't add all tets' faces in the face list
else
    [D] = pdist2(V,V(b,:),'euclidean');
end

if(max(max(D))==inf)
    error(['The vertices set is not connected!\n']);
end

if(false)
    RBF_D = 1./(log(D+10).*(D+10).*(D+10));
    norm_RBF_D = bsxfun(@times,RBF_D,1./sum(RBF_D,2));
    small_weight_region = (norm_RBF_D<0.005)*bc;

    handle_sample_dist_list_path = [folder_path,'\','dist_list.dmat'];
    handle_sample_dist_list = readDMAT(handle_sample_dist_list_path);
    threshold = 0.39;%handle_sample_dist_list(length(handle_sample_dist_list)-10)
    small_weight_region = (D>threshold)*bc;%280,330
else
    % nearest k handles as support
     D_rank_matrix = zeros(size(D));
     for i=1:1:size(D,1)
        [D_rank_matrix(i,:),~] = tiedrank(D(i,:));
     end
    LHS = (D_rank_matrix<=k);
end

