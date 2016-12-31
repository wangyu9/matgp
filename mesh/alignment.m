function [A] = alignment(source,dest,varargin)

% example usage:
% [Vg,Tg,Fg] = readMESH('unmapped.mesh');
% [Vc,Fc] = readOBJ('armadilloman.obj');
% A = alignment(Vg,Vc,'maxmin');
% V_aligned = [Vg,ones(size(Vg,1),1)]*A';
% V_aligned = V_aligned(:,1:3);

[method,bs,bd] = parse_optional_inputs(varargin);

switch(method)
    case 'maxmin'
        bs = maxmin_feature_point_indices(source);
        bd = maxmin_feature_point_indices(dest);
    case 'given'
        % use user input index
        assert(size(bs,1)>=3);
        assert(size(bd,1)>=3);
    otherwise
        error('Unknown alignment method!\n');
end

% assume there is a map from source to dest such that vd = A*vs.
% compute the pseudoinverse to get A

% we shoould have VD = VS * A', or VD'=A*VS'; if put all vs,vd vectors in rows
VD = dest(bd, :);
VD = [VD,ones(size(VD,1),1)];
VS = source(bs, :);
VS = [VS,ones(size(VS,1),1)];
% so (VD'*VS) = A*(VS'*VS)
A = (VD'*VS)/(VS'*VS);

end

function [D] = maxmin_feature_points(V)
    b = maxmin_feature_point_indices(V);
    D = V(b,:);
end

function [b] = maxmin_feature_point_indices(V)
    [~,bxmin] = min(V(:,1));
    [~,bxmax] = max(V(:,1));
    [~,bymin] = min(V(:,2));
    [~,bymax] = max(V(:,2));
    [~,bzmin] = min(V(:,3));
    [~,bzmax] = max(V(:,3));
    b = [bxmin;bxmax;bymin;bymax;bzmin;bzmax;];
end

function [method,bs,bd] = parse_optional_inputs(remaining_inputs)
  % parse remaining optional inputs
  ii = 1;
  % number of remaining inputs
  num_ri = numel(remaining_inputs);
  while(ii <= numel(remaining_inputs))
    switch remaining_inputs{ii}
%         case 'method'
%             ii = ii + 1;
%             method = remaining_inputs{ii};
%             switch method
%                 case ''
%             end
        case 'maxmin'
            method = 'maxmin';
        case 'given'
            method = 'given';
            assert(ii+2<=num_ri);
            bs = remaining_inputs{ii+1};
            bd = remaining_inputs{ii+2};
            ii = ii + 2;
        otherwise
            error('unknown option type!\n');            
    end
    ii = ii + 1;
  end
end