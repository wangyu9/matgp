%%
theta = [1/3,1/3,0]*pi;
theta(3) = pi - sum(theta);
%%
A = [0,cos(theta(2)),cos(theta(3));cos(theta(1)),0,cos(theta(3));cos(theta(1)),cos(theta(2)),0]
det(A)
%%
V = [0,0,0;1,0,0;-0.4,0.8,0];
F = [1,2,3];
B = eye(3);
%%
%VH = V;
%FH = F;
%bf = 1:size(F,1);
bf = find(V(F(:,1),1)<40);
%bf = 1:30;
VH = [V(F(bf,1),:);V(F(bf,2),:);V(F(bf,3),:);];
FH = (1:length(bf))'+length(bf)*ones(length(bf),1)*[0,1,2];
BH = eye(size(VH,1));
%BFH = eye(size(V,1));
%
for i=1:5
VH_old = VH;
FH_old = FH;
[VH,FH,MH] = upsample_with_faces_index(VH,FH,'KeepDuplicates',true);
%[~,~,MFH] = upsample_with_faces_index(VH_old,FH_old,'Interpolate',false);
BH = MH*BH;
%BFH = MFH*BFH;
end

%
%cl = [0,0,0;0,0,0;0,0,0];
%cq = [0,0,0;0,1,0;0,0,0];
%
[G] = grad(VH,FH);
fh = size(FH,1);
u = zeros(size(BH,1),1);
f = length(bf);
for i=1:length(bf)
    tid = bf(i);
    for j=1:3
        u = u + BH(:,(j-1)*f+i) * cl(tid,j);
        jp = mod(j-1-1,3)+1;
        jn = mod(j+1-1,3)+1;
        u = u + BH(:,(jp-1)*f+i) .* BH(:,(jn-1)*f+i) * cq(tid,j);
    end
end
%u = BH(:,1*3+0).*BH(:,2*3+0);
tid = 0;
bid = 1;
%u = BH(:,bid*3+tid);
u = u/(max(u)-min(u)) * max(max(VH)-min(VH));
u = u - min(u);
render_mesh3([VH(:,1),u,VH(:,2)],FH);%,'EdgeColor',[0,0,0]);
%hold on;
%render_mesh3([VH(:,1),zeros(size(VH,1),1),VH(:,2)],F,'EdgeColor',[0,0,0]);
%%
U = G*u;
U = reshape(U,[fh,3]);
u = U(:,2);
% 
render_mesh3(VH,FH,'EdgeColor',[0,0,0],'FaceColor',u);
% embree_render_mesh(V,F,)
%trimesh(F,V(:,1),V(:,2),V(:,3),u);
%tsurf(F,V(:,1:2),'EdgeColor',[0,0,0],'FaceColor',value2color((u-min(u)/(max(u)-min(u)))));
%%
%%
V = [0,0,0;1,0,0;-0.4,0.8,0;-0.4,-0.8,0];
F = [1,2,3;1,4,2];
%%
B = eye(3);
render_mesh3(V,F);
