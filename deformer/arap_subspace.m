function [Var,data,SS,R] = arap_subspace(varargin)
  % ARAP Solve for the as-rigid-as-possible deformation according to various
  % manifestations including:
  %   (1) "As-rigid-as-possible Surface Modeling" by [Sorkine and Alexa 2007]
  %   (2) "A local-global approach to mesh parameterization" by [Liu et al.
  %     2010] or "A simple geometric model for elastic deformation" by [Chao et
  %     al.  2010]
  %   (3) Adapted version of "As-rigid-as-possible Surface Modeling" by
  %     [Sorkine and Alexa 2007] presented in section 4.2 of or "A simple
  %     geometric model for elastic deformation" by [Chao et al.  2010]
  %
  % U = arap(V,F,b,bc) given a rest mesh (V,F) and list of constraint vertex
  % indices (b) and their new postions (bc) solve for pose mesh positions (U),
  % using default choice for 'Energy'
  %
  % U = arap(V,F,b,bc,'ParameterName','ParameterValue',...)
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by dim list of constraint positions for b
  %   Optional:
  %     'Energy'
  %       followed by a string specifying which arap energy definition to use.
  %       One of the following:
  %         'spokes'  "As-rigid-as-possible Surface Modeling" by [Sorkine and
  %           Alexa 2007], rotations defined at vertices affecting incident
  %           edges
  %         'elements'  "A local-global approach to mesh parameterization" by
  %           [Liu et al.  2010] or "A simple geometric model for elastic
  %           deformation" by [Chao et al.  2010], rotations defined at
  %           elements (triangles or tets) 
  %         'spokes-and-rims'  Adapted version of "As-rigid-as-possible Surface
  %           Modeling" by [Sorkine and Alexa 2007] presented in section 4.2 of
  %           or "A simple geometric model for elastic deformation" by [Chao et
  %           al.  2010], rotations defined at vertices affecting incident
  %           edges and opposite edges
  %     'Var0' #Var by dim list of initial guess variables
  %       dim by dim by #C list of linear transformations initial guesses,
  %       optional (default is to use identity transformations)
  %     'Data' see output
  %     'Groups'
  %       followed by #V list of group indices (1 to k) for each vertex, such 
  %       that vertex i is assigned to group G(i)
  %     'Tol'
  %       stopping critera parameter. If variables (linear transformation matrix
  %       entries) change by less than 'tol' the optimization terminates,
  %       default is 0.75 (weak tolerance)
  %     'MaxIter'
  %       max number of local-global iterations, default is 10
  %     'Dynamic' 
  %        #V by dim list of external forces
  %     'TimeStep'
  %        scalar time step value
  %     'Vm1' 
  %        #V by dim positions at time t-1
  %     'Tikhonov' followed by constant Tikhonov regularization parameter
  %       alpha:
  %       http://en.wikipedia.org/wiki/Tikhonov_regularization#Relation_to_probabilistic_formulation
  %     'Flat' followed by whether to add the constraint that Z=0 {false}
  %     'RemoveRigid' followed by whether to add a constraint that places an
  %       arbitrary point at the origin and another along the x-axis {false}
  % Outputs:
  %   U  #V by dim list of new positions
  %   data  struct of reusable data
  %     .CSM dim*n by dim*n sparse matrix containing special laplacians along the
  %       diagonal so that when multiplied by repmat(U,dim,1) gives covariance
  %       matrix elements, can be used to speed up next time this function is
  %       called, see function definitions
  %
  % Known issues: 'Flat',true + 'Energy','elements' should only need a 2D
  % rotation fit, but this does 3D to stay general (e.g. if one were to use
  % 'spokes' then the edge-set cannot be pre-mapped to a common plane)
  %
  % See also: takeo_arap
  %

  % parse input
  V = varargin{1};
  F = varargin{2};
  b = varargin{3};
  bc = varargin{4};
  G = [];

  % number of vertices
  n = size(V,1);
  % number of elements (tets or triangles)
  tf = size(F,1);
  
  assert(isempty(b) || max(b) <= n);
  assert(isempty(b) || min(b) >= 1);

  indices = 1:n;
  max_iterations = 100;
  tol = 0.001;
  interior = indices(~ismember(indices,b));% interior will be the complementary part of b, in 1:n
  Var = [];% U = [];
  Vm1 = [];
  % default is Sorkine and Alexa style local rigidity energy
  energy = 'spokes';
  % default is no external forces
  fext = [];
  % defaults is unit time step
  h = 1;
  % Tikhonov regularization alpha
  alpha_tik = 0;
  % flatten/parameterization
  flat = false;
  % remove rigid transformation invariance
  remove_rigid = false;
  G = [];
  debug = false;
  data = [];

  bases = [];
  cluster = [];
  
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Energy','Var0','Data','Tikhonov','Groups','Tol', ...
     'MaxIter','Dynamic','TimeStep','Vm1','Flat','RemoveRigid','Debug','Bases','Cluster'}, ...
    {'energy','Var0','data','alpha_tik','G','tol','max_iterations','fext', ...
    'h','Vm1','flat','remove_rigid','debug','bases','cluster'});
  v = 5;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if exist('Var0','var')
    Var = Var0;
  end
  
  if isempty(fext)
    dynamic = false;
    fext = zeros(size(V));
  else
    dynamic = true;
  end


  if isempty(bases)
     bases = speye(n,n); 
  end
  assert(size(bases,1)==n); % bases is of the dimensions of n*m
  m = size(bases,2); % m is the number of bases.
  
  if isempty(G)
    k = n;
  else
    k = max(G);
  end
  if flat
    [ref_V,ref_F,ref_map] = plane_project(V,F);
    assert(strcmp(energy,'elements'),'flat only makes sense with elements');
  else
    ref_map = 1;
    ref_V = V;% rest pose
    ref_F = F;
  end
  dim = size(ref_V,2);
  if isempty(bc)
    bc = sparse(0,dim);
  end

  assert(dim == size(bc,2));

  if isempty(Var)
    %if(dim == 2) 
    %  U = laplacian_mesh_editing(V,F,b,bc);
    %else
    %  U = V;    
    %end
    Var = zeros(m,dim); % initial guess
  end
  assert(m == size(Var,1));
  assert(dim == size(Var,2));

  if dynamic
    Var0 = Var(:,1:dim);
    if isempty(Vm1)
      Vm1 = Var0;
    end
    M = massmatrix(V,F);
    DQ = 0.5*1/h^2*bases'*M*bases;
    vel = (Var0-Vm1)/h;
    Dl = 1/(h^2)*bases'*M*bases*(-Var0 - h*vel) - fext;
    %Dl = 1/h^3*M*(-2*V0 + Vm1) - fext;
  else
    DQ = sparse(m,m);
    Dl = sparse(m,dim);
  end

  if ~isfield(data,'Q') || isempty(data.Q)
    data.Q = bases' * cotmatrix(V,F) * bases;
  end

  rr.b = cell(dim,1);
  rr.bc = cell(dim,1);
  if ~isfield(data,'rr')
    data.rr = [];
    if ~isfield(data.rr,'preF')
      data.rr.preF = cell(dim,1);
    end
  end
  if ~isfield(data,'preF')
    data.preF = [];
  end
  if remove_rigid
    error('The following has not been updated/implemented for subspace arap');
    if ~isempty(b)
      warning('RemoveRigid`s constraints are not typically wanted if |b|>0');
    end
    % the only danger is picking two points which end up mapped very close to
    % each other
    [~,f] = farthest_points(V,dim);
    for c = 1:dim
      rr.b{c} = f([1:dim-(c-1)])';
      rr.bc{c} = zeros(numel(rr.b{c}),1);
    end
    % Immediately remove rigid transformation from initial guess
    U = bsxfun(@minus,U,U(f(1),:));
    % We know that f(1) is at origin
    switch dim
    case 3
      % rotate about y-axis so that f(2) has (x=0)
      theta = atan2(U(f(2),1),U(f(2),3));
      R = axisangle2matrix([0 1 0],theta);
      U = U*R;
      % rotate about x-axis so that f(2) has (y=0)
      theta = atan2(U(f(2),3),U(f(2),2))+pi/2;
      R = axisangle2matrix([1 0 0],theta);
      U = U*R;
      % rotate about z-axis so that f(3) has (x=0)
      theta = atan2(U(f(3),2),U(f(3),1))+pi/2;
      R = axisangle2matrix([0 0 1],theta);
      U = U*R;
    case 2
      % rotate so that f(2) is on y-axis (x=0)
      theta = atan2(U(f(2),2),U(f(2),1))+pi/2;
      R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
      U = U*R;
    otherwise
      error('Unsupported dimension');
    end
  end

  all = [interior b];

  % cholesky factorization
  %cholL = chol(-L(interior,interior),'lower');
  % Why lu and not cholesky?
  %[luL,luU,luP,luQ,luR] = lu(-L(interior,interior));

  %R = repmat(eye(dim,dim),[1 1 n]);

  if ~isfield(data,'ae') || isempty(data.ae)
    data.ae = avgedge(V,F);
  end

  % build covariance scatter matrix used to build covariance matrices we'll
  % later fit rotations to
  if ~isfield(data,'CSM') || isempty(data.CSM)
    assert(size(ref_V,2) == dim);
    data.CSM = covariance_scatter_matrix(ref_V,ref_F,'Energy',energy) * repdiag(bases,dim);
    if flat
      data.CSM = data.CSM * repdiag(ref_map',dim);
    end
    
    % if there are groups then condense scatter matrix to only build
    % covariance matrices for each group
    if ~isempty(G)
      if strcmp(energy,'elements') && numel(G) ~= size(F,1) && numel(G) == n
        % groups are defined per vertex, convert to per face using mode
        G = mode(G(F),2);
      end
      G_sum = group_sum_matrix(G,k);
      %CSM = [G_sum sparse(k,n); sparse(k,n) G_sum] * CSM;
      data.CSM = repdiag(G_sum,dim) * data.CSM;
    end
  end

  % precompute rhs premultiplier
  if ~isfield(data,'K') || isempty(data.K)
    [~,data.K] = arap_rhs(ref_V,ref_F,[],'Energy',energy);
    data.K = repdiag(bases',dim) * data.K;
    if flat
      error('The following has not been updated!');
      data.K = repdiag(ref_map,dim) * data.K;
    end
  end

  % initialize rotations with identies (not necessary)
  R = repmat(eye(dim,dim),[1 1 size(data.CSM,1)/dim]);

  iteration = 0;
  Var_prev = Var;
  data.energy = inf;
  while true

    if iteration > max_iterations
      if debug
        fprintf('arap: Iter (%d) > max_iterations (%d)\n',iteration,max_iterations);
      end
      break;
    end

    if iteration > 0
      change = max(abs(Var(:)-Var_prev(:)));
      if debug
        fprintf('arap: iter: %d, change: %g, energy: %g\n', ...
          iteration,change,data.energy);
      end
      if change <tol*data.ae
        if debug
          fprintf('arap: change (%g) < tol*ae (%g * %g)\n',change,tol,data.ae);
        end
        break;
      end
    end

    Var_prev = Var;

    % energy after last global step
    Var(b,:) = bc;
    %E = arap_energy(V,F,U,R);
    %[1 E]

    % compute covariance matrix elements
    S = zeros(size(data.CSM,1),dim);
    S(:,1:dim) = data.CSM*repmat(Var,dim,1);
    % dim by dim by n list of covariance matrices
    SS = permute(reshape(S,[size(data.CSM,1)/dim dim dim]),[2 3 1]);
    % fit rotations to each deformed vertex
    R = fit_rotations(SS,'SinglePrecision',false);

    % for debug purpose
%    R = repmat(eye(dim,dim),[1 1 size(data.CSM,1)/dim]);
    
    % energy after last local step
    Var(b,:) = bc;
    %E = arap_energy(V,F,U,R);
    %[2 E]


    % This is still the SLOW way of building the right hand side, the entire
    % step of building the right hand side should collapse into a
    % #handles*dim*dim+1 by #groups matrix times the rotations at for group
    
    % distribute group rotations to vertices in each group
    if ~isempty(G)
      R = R(:,:,G);
    end

    Var(b,:) = bc;

    %B = arap_rhs(V,F,R);
    Rcol = reshape(permute(R,[3 1 2]),size(data.K,2),1);
    Bcol = data.K * Rcol;
    B = reshape(Bcol,[size(Bcol,1)/dim dim]);

    % RemoveRigid requires to solve each indepently
    if remove_rigid
      error('Not updated this part!');
      for c = 1:dim
        eff_b = [b rr.b{c}];
        eff_bc = [bc(:,c);rr.bc{c}];
        [U(:,c),data.rr.preF{c}] = min_quad_with_fixed( ...
          -0.5*data.L+DQ+alpha_tik*speye(size(data.L)), ...
          -B(:,c)+Dl(:,c),eff_b,eff_bc,[],[],data.rr.preF{c});
      end
    else 
      [Var,data.preF] = min_quad_with_fixed( ...
        -0.5*data.Q+DQ+alpha_tik*speye(size(data.Q)), ...
        -B+Dl,b,bc,[],[],data.preF);
    end
    energy_prev = data.energy;
    data.energy = trace(Var'*(-0.5*data.Q)*Var+Var'*(-B)+Var0'*(-0.5*data.Q*Var0));
    if data.energy > energy_prev
      %if debug
        fprintf('arap: energy (%g) increasing (over %g)\n', ...
          data.energy,energy_prev);
      %end
      break;
    end
    %U(interior,:) = -L(interior,interior) \ (B(interior,:) + L(interior,b)*bc);
    %U(interior,:) = -L(interior,all)*L(all,interior) \ (L(interior,all)*B(all,:) + L(interior,all)*L(all,b)*bc);
    %U(interior,:) = luQ*(luU\(luL\(luP*(luR\(B(interior,:)+L(interior,b)*bc)))));
    %U(interior,:)=cholL\((B(interior,:)+L(interior,b)*bc)'/cholL)';
    iteration = iteration + 1;
  end
  Var = Var(:,1:dim);
end
