%% Generate Some Synthetic Data
clear all

nInstances = 200;%200
nVars = 250;%250
sparsityFactor = .5;
flipFactor = .1;
X = [ones(nInstances,1) randn(nInstances,nVars-1)];
w = randn(nVars,1).*(rand(nVars,1) < sparsityFactor);
y = sign(X*w);
flipPos = rand(nInstances,1) < flipFactor;
y(flipPos) = -y(flipPos);
        
%% Set up optimization problem
w_init = zeros(nVars,1);

lambda = 1;
lambdaVect = lambda*[0;ones(nVars-1,1)];

funObj = @(w)LogisticLoss(w,X,y);

%% Set Optimization Options
gOptions.maxIter = 2000;
gOptions.verbose = 1; % Set to 0 to turn off output
options.corrections = 10; % Number of corrections to store for L-BFGS methods

%% Run Solvers

fprintf('Spectral Projected Gradient\n');
options = gOptions;
wSPG = L1General2_SPG(funObj,w_init,lambdaVect,options);
pause;

fprintf('Barzilai-Borwein Soft-Threshold\n');
options = gOptions;
wBBST = L1General2_BBST(funObj,w_init,lambdaVect,options);
pause;

fprintf('Barzilai-Borwein Sub-Gradient\n');
options = gOptions;
wBBSG = L1General2_BBSG(funObj,w_init,lambdaVect,options);
pause;

fprintf('Optimal Projected Gradient\n');
options = gOptions;
options.L = 1/size(X,1);
wSPG = L1General2_OPG(funObj,w_init,lambdaVect,options);
pause;

fprintf('Diagonally-Scaled Soft-Threshold\n');
options = gOptions;
funObj_DiagHess = @(w)LogisticLossDiagHess(w,X,y);
wDSST = L1General2_DSST(funObj_DiagHess,w_init,lambdaVect,options);
pause;

fprintf('Orthant-Wise Learning\n');
options = gOptions;
options.quadraticInit = 1;
wOWL = L1General2_OWL(funObj,w_init,lambdaVect,options);
pause;

fprintf('Active Set\n');
options = gOptions;
wAS = L1General2_AS(funObj,w_init,lambdaVect,options);
pause;

fprintf('Two-Metric Projection\n');
options = gOptions;
wTMP = L1General2_TMP(funObj,w_init,lambdaVect,options);
pause;

fprintf('Projected Scaled Sub-Gradient (Gafni-Bertsekas variant)\n');
options = gOptions;
wPSSgb = L1General2_PSSgb(funObj,w_init,lambdaVect,options);
pause;

fprintf('Projected Scaled Sub-Gradient (Sign-Projection variant)\n');
options = gOptions;
wPSSsp = L1General2_PSSsp(funObj,w_init,lambdaVect,options);
pause;

fprintf('Projected Scaled Sub-Gradient (Active-Set variant)\n');
options = gOptions;
wPSSas = L1General2_PSSas(funObj,w_init,lambdaVect,options);
