function [B] = farthest_point_sampling(V,T,F,b,num,seed,varargin)
% typical usage:
% 3D: [B] = farthest_point_sampling(V,T,F,[],60,[],'noboundary');% note the
% F must be there
% 2D: [B] = farthest_point_sampling(V,[],F,[],60,[],'noboundary');

% V,T,F is the vertices, tets and faces of mesh.
% b is the exising list of sampled vertices.
% num is the number of vertices to sample.
% seed is the index of vertice to start if b is empty.
if(size(b,1)==1)
    b = b';
end
assert(isempty(b)||size(b,2)==1);
assert(isempty(seed)||length(seed)==1);

% varibles from optional inputs
noboundary = false;
% parse optional inputs
ii = 1;
num_ri = numel(varargin);
while(ii <= numel(varargin))
   switch varargin{ii}
       case 'noboundary'
           noboundary = true;
       otherwise
           error(['''' varargin{ii} ''' is not a valid parameter']);
   end
   ii = ii + 1;
end
% end of parsing

if(noboundary)
    if(size(T,1)~=0)
        [boundary_list,~] = boundary_vertices(T);
        not_boundary_list = ones( max(max(T)), 1);
        not_boundary_list(boundary_list) = 0;
        not_boundary_list = find(not_boundary_list);
    else
        [boundary_list,~] = boundary_vertices(F);
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

new_D = [];
%new_one = [];
D = [];
D_closet = Inf(n,1); % a vector storing the distance to the closet point in B.

if(noboundary)
    D_closet(boundary_list) = 0;
end



% this is a do-while loop in c. Matlab is so stupid to not have it.
i = 1;
flag = true;
while flag
   % do stuff
   
   B = [B;new_B];
   [new_D] = mesh_dist(V,T,F,new_B);
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

assert(length(B)==length(unique(B)));

end

function [D] = mesh_dist(V,T,F,b)
    if(false)% TODO
        [D] = geo_dist(V,T,F,b);% add T here since tetgen won't add all tets' faces in the face list
    else
        [D] = pdist2(V,V(b,:),'euclidean');
    end
end