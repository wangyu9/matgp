function [V,F] = quiver_mesh(normal,scale)%start,normal,scale,varargin)

% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
% GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;...
%               norm(cross(A,B)) dot(A,B)  0;...
%               0              0           1];
% 
% FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

% UU = @(Fi,G) Fi*G*inv(Fi);


a=[1 0 0]'; b=normal';
% U = UU(FFi(a,b), GG(a,b))
U = fcn_RotationFromTwoVectors(a, b);

m = 1;
n = 10;
[W,G] = rectangle_mesh(-1,-1,2/m,2/n,m,n);
R = 0.1;

V = [W(:,1),R*sin(pi*W(:,2)),R*cos(pi*W(:,2))];
F = G;


m = 1;
n = 10;
[W,G] = rectangle_mesh(-1,-1,2/m,2/n,m,n);


V2 = [(W(:,1)+1)*0.3+1,R*(1-W(:,1)).*sin(pi*W(:,2)),R*(1-W(:,1)).*cos(pi*W(:,2))];
F2 = G;

[V,~,F] = merge_mesh(V,[],F,V2,[],F2);

V = (V-[-1,0,0])*U'*scale;

function R=fcn_RotationFromTwoVectors(v1, v2)
% R*v1=v2
% v1 and v2 should be column vectors and 3x1

% 1. rotation vector
w=cross(v1,v2);

if norm(w)<1e-15
    R = eye(size(v1,1));
else
    w=w/norm(w);
    w_hat=fcn_GetSkew(w);
    % 2. rotation angle
    cos_tht=v1'*v2/norm(v1)/norm(v2);
    tht=acos(cos_tht);
    % 3. rotation matrix, using Rodrigues' formula
    R=eye(size(v1,1))+w_hat*sin(tht)+w_hat^2*(1-cos(tht));
end

function x_skew=fcn_GetSkew(x)
x_skew=[0 -x(3) x(2);
 x(3) 0 -x(1);
 -x(2) x(1) 0];

