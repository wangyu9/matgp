function W = biharmonic_bounded(varargin)
  % BIHARMONIC_BOUNDED Compute biharmonic bounded coordinates, using quadratic
  % optimizer
  %
  % W = biharmonic_bounded(V,F,b,bc)
  % W = biharmonic_bounded(V,F,b,bc,'ParameterName',ParameterValue)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %    Optional:
  %      'k'  followed by integer, solve k-harmonic problem {2}
  %      'QuadProgParam'  followed by mosek param struct, see optimset. If not
  %        given then *some* default options are used
  %      'OptType'  of optimizer to use {best available}:
  %        'quad'
  %        'least-squares'
  %        'conic'
  %      'POU'  true or false, enforce partition of unity explicitly {false}
  %      'Low'  lower bound {0}
  %      'Up'  upper bound {1}
  %      'ShapePreserving'  #V rigidity mask, where 0 means not enforced, 1 means strongly
  %        enforced
  %  
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions
  %

  V = varargin{1};
  F = varargin{2};
  b = varargin{3};
  bc = varargin{4};

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % default options
  pou = false;
  k = 2;
  % check for mosek and set its parameters
  [param,mosek_exists] = default_quadprog_param();
  param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
  if mosek_exists
    opt_type = 'conic';
  else
    opt_type = 'quad';
  end
  % default bounds
  low = 0;
  up = 1;
  R = sparse(n,1);

  % parse optional inputs
  ii = 5;
  while(ii <= nargin)
    switch varargin{ii}
    case 'POU'
      if (ii+1)<=nargin && islogical(varargin{ii+1})
        ii = ii + 1;
        assert(ii<=nargin);
        pou = varargin{ii};
      else
        pou = true;
      end
    case 'k'
      ii = ii + 1;
      assert(ii<=nargin);
      k = varargin{ii};
      if k > 2
          opt_type = 'quad';
      end
    case 'QuadProgParam'
      ii = ii + 1;
      assert(ii<=nargin);
      param = varargin{ii};
    case 'Low'
      ii = ii + 1;
      assert(ii<=nargin);
      low = varargin{ii};
    case 'Up'
      ii = ii + 1;
      assert(ii<=nargin);
      up = varargin{ii};
    case 'OptType'
      ii = ii + 1;
      assert(ii<=nargin);
      opt_type = varargin{ii};
    case 'ShapePreserving'
      ii = ii + 1;
      assert(ii<=nargin);
      R = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii+1;
  end

  % Build discrete laplacian and mass matrices used by all handles' solves
  if(size(F,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,F);
    M = massmatrix3(V,F,'barycentric');
  else
    %L = cotmatrix(V,F);
    %M = massmatrix(V,F,'voronoi');
    L = cotmatrix_embedded(V,F);
    M = massmatrix_embedded(V,F,'voronoi');
  end
  % NORMALIZE MASSMATRIX (THIS IS IMPORTANT!!)
  M = M./max(abs(diag(M)));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP SOLVER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if any(R(:)) && ~strcmp(opt_type,'quad')
    error( [ ...
      'Enforcing shape preserving regions only supported in' ...
      ' conjunction with opt_type=''quad''']);
  end

  %% Shape preserving laplacian
  %Lsp = bsxfun(@times,L,sparse(R)*shape_preserving_weight);
  %spy(Lsp)
  shape_preserving_weight = 100;
  E = edges(F);
  ne = size(E,1);
  % Build gradient
  G = sparse( ...
    [1:ne 1:ne]', ...
    [E(:,1);E(:,2)], ...
    shape_preserving_weight*[mean(R(E),2);-mean(R(E),2)], ...
    ne, ...
    n);
  Lsp = G' * G;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP PROBLEM AND SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(pou)
    % Enforce partition of unity as explicity constraints: solve for weights
    % of all handles simultaneously
    if(strcmp(opt_type,'quad'))
      % biharmonic system matrix
      Qi = -L;
      for ii = 2:k
        Qi = -L*(M\Qi);
      end
      if k>3
        Qi = (Qi + Qi')*0.5;
      end
      Qi = Qi + Lsp;
      Q = sparse(m*n,m*n);
      % Q is sparse matrix with Qi along diagonal
      for ii = 1:m
        d = (ii - 1)*n + 1;
        Q(d:(d + n-1), d:(d + n-1)) = Qi;
      end
      % linear constraints: partition of unity constraints and boundary
      % conditions
      PA = repmat(speye(n,n),1,m);
      Pb = ones(n,1);
      % boundary conditions
      BCAi = speye(n,n);
      BCAi = BCAi(b,:);
      BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
      % BCA is sparse matrix with BCAi along diagonal
      for ii = 1:m
        di = (ii - 1)*size(BCAi,1) + 1;
        dj = (ii - 1)*size(BCAi,2) + 1;
        BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
      end
      BCb = bc(:);
      % set bounds
      ux = up.*ones(m*n,1);
      lx = low.*ones(m*n,1);
      if(mosek_exists)
        fprintf('Quadratic optimization using mosek...\n');
      else
        fprintf('Quadratic optimization using matlab...\n');
      end
      fprintf( [ ...
        '  minimize:     x''LM\\Lx\n' ...
        'subject to: %g <= x <= %g, ???_i xi = 1\n'], ...
        low,up);
      tic;
      W = quadprog(Q,zeros(n*m,1),[],[],[PA;BCA],[Pb;BCb],lx,ux,[],param);
      toc
      W = reshape(W,n,m);
    else
      error( [ ...
        'Enforcing partition of unity only support in conjunction with ' ...
        'opt_type=''quad''']);
    end
  else
    % Drop partition of unity constraints, solve for weights of each handle
    % independently then normalize to enforce partition of unity
    if(strcmp(opt_type,'quad'))
      % build quadratic coefficient matrix (bilaplacian operator)
      Q = -L;
      for ii = 2:k
        Q = -L*(M\Q);
      end
      if k>=3
        Q = (Q + Q')*0.5;
      end
      % set bounds
      ux = up.*ones(n,1);
      lx = low.*ones(n,1);
    elseif(strcmp(opt_type,'least-squares'))
      % solve same problem but as least-squares problem see mosek documention
      % for details
      I = speye(n);
      Z = sparse(n,n);
      Q = [Z,Z;Z,I];
      assert(k == 2);
      F = sqrt(M)\L;
      c = zeros(n,1);
      B = [F,-I];
      ux = [up.*ones(n,1) ;  Inf*ones(n,1)];
      lx = [low.*ones(n,1); -Inf*ones(n,1)];
    elseif(strcmp(opt_type,'conic'))
      % solve same problem but as conic problem see mosek documention for
      % details
      %
      % Variable names are consistent with mosek doc 7.9.1:
      % [x z t]
      assert(k==2);
      F = sqrt(M)\L;
      prob.c = [zeros(2*n,1); 1];
      I = speye(n);
      prob.a = [F,-I,zeros(n,1)];
      prob.blc = zeros(n,1);
      prob.buc = zeros(n,1);
      prob.bux = [ up.*ones(n,1);  Inf*ones(n,1);  Inf];
      prob.blx = [ low.*ones(n,1); -Inf*ones(n,1); 0];
      prob.cones = cell(1,1);
      prob.cones{1}.type = 'MSK_CT_QUAD';
      t_index = 2*n +1;
      z_indices = (n+1):(2*n);
      prob.cones{1}.sub = [t_index z_indices];
    else
      error('Bad opt_type');
    end

    % number of handles
    m = size(bc,2);
    % allocate space for weights
    W = zeros(n,m);
    tic;
    % loop over handles
    for i = 1:m
      if(strcmp(opt_type,'quad'))
        % enforce boundary conditions via lower and upper bounds
        %lx(b) = bc(:,i);
        %ux(b) = bc(:,i);
        Aeq = speye(n,n);
        Aeq = Aeq(b,:);
        if(mosek_exists)
          fprintf('Quadratic optimization using mosek...\n');
        else
          fprintf('Quadratic optimization using matlab...\n');
        end
        fprintf( [ ...
          '  minimize:     x''LM\\Lx\n' ...
          'subject to: %g <= x <= %g\n' ], ...
          low,up);
        % if mosek is not available, then matlab will complain that sparse
        % matrices are not yet supported...
        [x,fval,err] = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,i),lx,ux,[],param);
        if(err ~= 1)
          fprintf([...
            '----------------------------------------------------------\n' ...
            'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
            ' solution may be inaccurate...\n' ...
            '----------------------------------------------------------\n' ...
            ]);
        end
        
        if false % this is removed by wangyu

        prob = [];
        [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(Q));
        prob.c = zeros(n,1);
        prob.a = sparse(0,n);
        prob.blc = sparse(0,1);
        prob.buc = sparse(0,1);
        prob.blx = lx;
        prob.bux = ux;
        % Enforce fixed values
        prob.bux(b) = bc(:,i);
        prob.blx(b) = bc(:,i);
        quiet = '';
        fprintf('Quadratic optimization using mosek...\n');
        [r,res]=mosekopt(['minimize' quiet],prob,param);
        % report_mosek_error(r,res);% this is removed by wangyu
        
        end

      elseif(strcmp(opt_type,'least-squares'))
        % enforce boundary conditions via lower and upper bounds
        lx(b) = bc(:,i);
        ux(b) = bc(:,i);
        fprintf('Quadratic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:       z''z\n' ...
          '  subject to: M\\Lx - z = 0\n' ...
          '  and          %g <= x <= %g\n'], ...
          low,up);
        x = quadprog(Q,zeros(2*n,1),[],[],B,c,lx,ux,[],param);
      elseif(strcmp(opt_type,'conic'))
        prob.bux(b) = bc(:,i);
        prob.blx(b) = bc(:,i);
        fprintf('Conic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:         t\n' ...
          '  subject to: M\\Lx - z = 0,\n' ...
          '             t >= sqrt(z''z),\n' ...
          '               %f <= x <= %f\n'], ...
          low,up);
        [r,res]=mosekopt('minimize echo(0)',prob,param);
        % check for mosek error
        if(r == 4006)
          warning(['MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg ...
            'The solution is probably OK, but ' ...
            'to make this error go away, increase: ' ...
            'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
            n]);
        elseif(r ~= 0)
          error(['FATAL MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg]);
        end
        % extract solution from result
        x = res.sol.itr.xx;
      end
      % set weights to solution in weight matrix
      W(:,i) = x(1:n);
      fprintf('Lap time: %gs\n',toc);
    end
    t = toc;
    fprintf('Total elapsed time: %gs\n',t);
    fprintf('Average time per handle: %gs\n',t/m);
  end

end
