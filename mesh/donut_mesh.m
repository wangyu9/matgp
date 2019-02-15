function [V,F] = donut_mesh(VL,m,d,close_loop)
%%
% See example_donut_mesh.m
%%
%
% close_loop = false;
%%
n = size(VL,1);
VL = bsxfun(@plus,VL,[d,0,0]);
%%
%m = 80;
V = zeros(n*m,3);
for i= 1:m
    theta = 2*pi/m*i;
    R = blkdiag([cos(theta),-sin(theta);sin(theta),cos(theta)],1);
    V((1:n)+(i-1)*n,:) = VL*R';
end
%%
scatter3(V(:,1),V(:,2),V(:,3));
%%
N = m*n;
%F = zeros(2*m*n,3);
F = [];
if ~close_loop
    for i=1:m  
        Ia = (1:n-1)'+(i-1)*n;
        Ib = [(2:n)]'+(i-1)*n;
        Ip = mod(Ia+n-1,N)+1; % since matlab is 1-indexing...
        Iq = mod(Ib+n-1,N)+1;
        F = [F;[Ia,Ip,Ib];[Ib,Ip,Iq]];
    end    
else
    for i=1:m
    %for i=1:m-1% cut
    %for i=1:m/2-1
    %for i=m/2+1:m-1    
        Ia = (1:n)'+(i-1)*n;
        Ib = [(2:n),1]'+(i-1)*n;
        Ip = mod(Ia+n-1,N)+1; % since matlab is 1-indexing...
        Iq = mod(Ib+n-1,N)+1;
        F = [F;[Ia,Ip,Ib];[Ib,Ip,Iq]];
    end
end

