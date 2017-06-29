function [D] = earth_mover_distance(X,T,b)

%% Compute LB eigenstuff

nEigs = 100;
FEM = firstOrderFEM(X,T); % set up FEM matrices
[evecs,evals] = eigs(FEM.laplacian,FEM.vtxInnerProds,nEigs,'sm'); %original test was 300
evals = diag(evals);

%% Precompute quantities that can be useful for multiple EMD computations

structure = precomputeEarthMoversADMM(X, T, evecs(:,2:end));

%% Design distributions as centered around points for testing

%point0 = 1779; %nose
%point1 = 7943; %tail
point0 = b;

% Construct extrinsic Gaussians around these two points
h = mean( max(X) - min(X));
a0 = 2*(3.7651*h)^2;%50; % determines standard deviation
a1 = a0;%50;

dist0 = sqrt(sum(bsxfun(@minus,X,X(point0,:)).^2,2));
rho0 = exp(-dist0.^2*a0);
rho0 = rho0/sum(rho0);
showDescriptor(X,T,rho0); % displays mesh with a function on vertices
title('Rho 0');

%%
n = size(X,1);
D = zeros(n,1);

for i = 1:n

    dist1 = sqrt(sum(bsxfun(@minus,X,X(i,:)).^2,2));
    rho1 = exp(-dist1.^2*a1);
    rho1 = rho1/sum(rho1);
    %showDescriptor(X,T,rho1);
    %title('Rho 1');

    % Compute EMD

    [distance,J] = earthMoversADMM(X, T, rho0, rho1, structure);
    % fprintf('EMD = %g\n',distance);

    D(i) = distance;
    % Display momentum field

    % showFaceVectorField(X, T, J{1});
    % title('Momentum vector field');

end