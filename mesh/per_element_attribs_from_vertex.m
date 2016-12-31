function [A_element] = per_element_attribs_from_vertex(A,TF,type)
    % A is per vertex attributes
    % A_element is per element(trangles/tets)
    
    n = size(A,1);
    tf = size(TF,1);
    dc = size(TF,2);
    
    A_element = zeros(tf,size(A,2));

    switch(type)
        case 'uniform'
        % case 1: uniform average
            for i=1:1:tf
                for c=1:1:dc
                    A_element(i,:) = A_element(i,:)+A(TF(i,c),:);
                end
            end
            A_element = A_element / dc;

        % case 2: maybe cotan weighting, etc
            % TODO
            
        otherwise
            error(['Unknown type.\n']);
    end
end