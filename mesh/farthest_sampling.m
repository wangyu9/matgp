function [B] = farthest_sampling(D,num,varargin)


% varibles from optional inputs
b = [];
seed = 1;
include_seed = false;
% parse optional inputs
ii = 1;
num_ri = numel(varargin);
while(ii <= numel(varargin))
   switch varargin{ii}
       case 'b'
           b = varargin{ii+1};
           ii = ii + 1;
       case 'seed'
           seed = varargin{ii+1};
           ii = ii + 1;           
       otherwise
           error(['''' varargin{ii} ''' is not a valid parameter']);
   end
   ii = ii + 1;
end
% end of parsing

if(size(b,1)==1)
    b = b';
end
assert(isempty(b)||size(b,2)==1);
assert(isempty(seed)||length(seed)==1);


new_B = b;
remove_top = false;

if(length(b)==0)
    b = seed;
    num = num - 1;
    new_B = seed;
    % note: no need to update B here, will do it in the loop
end

B = [];

n = size(D,1);

new_D = [];
D_closet = Inf(n,1); % a vector storing the distance to the closet point in B.

% this is a do-while loop in c. Matlab is so stupid to not have it.
i = 1;
flag = true;
while flag
   % do stuff
   
   B = [B;new_B];
   [new_D] = D(:,B);
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
