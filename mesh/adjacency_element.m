function [Fb] = adjacency_element(F,b,threshold)

% find all elements such that at least have connections to threshold
% vertices in the list of b


IH = b;

    MAP = zeros(max(max(F)),1);
    MAP(IH,:) = 1;
    
    if(size(F,2)==3)
        numInIH = sum([MAP(F(:,1)),MAP(F(:,2)),MAP(F(:,3))],2);
    else
        assert(size(F,2)==4);
        numInIH = sum([MAP(F(:,1)),MAP(F(:,2)),MAP(F(:,3)),MAP(F(:,4))],2);
    end
    
    
    Fb = find(numInIH>=threshold);