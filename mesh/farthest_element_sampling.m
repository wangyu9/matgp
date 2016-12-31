function [B] = farthest_element_sampling(V,TF,b,num,seed,varargin)
% example
    %%if want the list of vertices of sampled elements, could do this:
    % bb = TF(B,:)';
    % bb = bb(:);
    % H = V(bb,:);
    % writeDMAT('H',H);    
    %%if want the group indexing further, could do this:
    % group = (1:1:size(TF,1));
    % group = repmat(group,size(TF,2),1);
    % group = group(:);
    % group = group -1; % make 0-index before save
    % writeDMAT('group',group);
    
    %%if want to know the distance to the closeast sampled element
% input
    % V,T,F is the vertices, tets and faces of mesh.
    % b is the exising list of sampled vertices.
    % num is the number of vertices to sample.
    % seed is the index of vertice to start if b is empty.
% output:
    % B is a list of sampled tet indices
    % IB is a list of sampled vertice indices  
    % I is a list of sampled 
    
    assert(isempty(b)||size(b,2)==1);
    assert(isempty(seed)||length(seed)==1);

    % varibles from optional inputs
    noboundary = false;
    dist_type = 'euclidean';%'geodesic';
    % parse optional inputs
    ii = 1;
    num_ri = numel(varargin);
    while(ii <= numel(varargin))
       switch varargin{ii}
           case 'noboundary'
               noboundary = true;
           case 'euclidean'
               dist_type = 'euclidean';
           case 'geodesic'
               dist_type = 'geodesic';
           otherwise
               error(['''' varargin{ii} ''' is not a valid parameter']);
       end
       ii = ii + 1;
    end
    % end of parsing

    if(noboundary)
        if(size(T,1)~=0)
            [~,boundary_list] = boundary_faces(T);% here is different from farthest_point_sampling
            not_boundary_list = ones( max(max(T)), 1);
            not_boundary_list(boundary_list) = 0;
            not_boundary_list = find(not_boundary_list);
        else
            [~,boundary_list] = boundary_edges(F);% here is different from farthest_point_sampling
            not_boundary_list = ones( max(max(F)), 1);
            not_boundary_list(boundary_list) = 0;
            not_boundary_list = find(not_boundary_list);
        end
    end

    if(exist('seed')~=1||length(seed)==0)
        if(noboundary)
            seed = not_boundary_list(1);
        else
            seed = 1;
        end
    end

    new_B = b;

    if(size(b,1)==0)
        b = seed;
        num = num - 1;
        new_B = seed;
        % note: no need to update B here, will do it in the loop
    end

    B = [];

    n = size(V,1);

    % here is different from farthest_point_sampling
    tf = size(TF,1);
    V_e = element_centers(V,TF);
%     if(size(T,1)~=0)
%         tf = size(T,1);
%         V_e = tet_centers(V,T);
%     else
%         tf = size(F,1);
%         V_e = face_centers(V,F);
%     end

    new_D = [];
    %new_one = [];
    D = [];
    D_closet = Inf(tf,1); % a vector storing the distance to the closet point in B.
    % here is different from farthest_point_sampling

    if(noboundary)
        D_closet(boundary_list) = 0;
    end



    % this is a do-while loop in c. Matlab is so stupid to not have it.
    i = 1;
    flag = true;
    while flag
       % do stuff

       B = [B;new_B];
       
       switch dist_type 
           case 'geodesic'
               [new_D] = tet_geo_dist(V,TF,new_B);
               error(['Seems to have some bugs.']);
           case 'euclidean'
               [new_D] = tet_dist(V_e,new_B);
       end     
       
       if(size(new_D,2)>1)
           new_D = min(new_D,[],2);% this together means min per row
       end
       update = find(new_D<D_closet);
       D_closet(update,:) = new_D(update,:); 

       [min_dist,new_B] = max(D_closet);
       if(min_dist==0)
           warning(['No possible sampling, exit.\n'])
           break;
       end

       % then it is part of do-while loop
       if(i>num)
        flag = false;
       end

       i = i + 1;
    end

    
    
end

function [D] = tet_dist(V,b)
   [D] = pdist2(V,V(b,:),'euclidean');
end

function [D] = tet_geo_dist(V,TF,b)
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

function [V_center] = element_centers(V,TF)
    V_center = per_element_attribs_from_vertex(V,TF,'uniform');
end

function [V_center] = tet_centers(V,T)
    V_center = element_centers(V,T);
end

function [V_center] = face_centers(V,F)
    V_center = element_centers(V,F);
end