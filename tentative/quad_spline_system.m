function [r] = quad_spline_system(V,F)
%%
%V = [0,0,0;1,0,0;-0.4,0.8,0;-0.4,-0.8,0];
%F = [1,2,3;1,4,2];
if false
%%
%V = [0,0,0;1,0,0;0.5,0.6,0;-0.2,-0.76,0;1.2,1.2,0];
%F = [1,2,3;1,4,2;2,5,3];
V = [0,0,0;1,0,0;0.5,0.6,0;-0.2,-0.76,0;1.2,1.2,0;2,0,0;1.3,-1,0];
F = [1,2,3;1,4,2;2,5,3;2,6,5;2,7,6;2,4,7];
%%
%[V,F] = readOBJ('woody.obj');
%V = V / max(abs(V));
[V,F] = readOBJ('alligator.obj');

%F = F(1:20,:);
end
%%
f = size(F,1);

dblA = doublearea(V,F);
EL = edge_lengths(V,F); % f x 3
EH = bsxfun(@times,dblA,1./EL);
DE1 = [F(:,2),F(:,3)];
DE2 = [F(:,3),F(:,1)];
DE3 = [F(:,1),F(:,2)];
DE = [DE1;DE2;DE3]; % 3f x 2
iF_DE = [(1:f)';(1:f)';(1:f)']; % index of which triangle each edge comes from.
iiF_DE = [ones(f,1);2*ones(f,1);3*ones(f,1)]; % in-triangle index of each edge
%
v1 = V(DE1(:,2),:) - V(DE1(:,1),:);
v2 = V(DE2(:,2),:) - V(DE2(:,1),:);
v3 = V(DE3(:,2),:) - V(DE3(:,1),:);
cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
Ecot = [cot23,cot31,cot12];
N = cross(v1,v2);
assert(min(N(:,3))>0);
%
cos12 = -dot(v1,v2,2)./(normrow(v1).*normrow(v2)); 
cos23 = -dot(v2,v3,2)./(normrow(v2).*normrow(v3)); 
cos31 = -dot(v3,v1,2)./(normrow(v3).*normrow(v1)); 
Ecos = [cos23,cos31,cos12];

% find boundary edges and interior edges. 
% build directed adjacency matrix
[~,iIDE,~] = intersect(DE,[DE(:,2),DE(:,1)],'rows');
bIDE = setdiff(1:size(DE,1),iIDE);
iIDE = iIDE(:);
bIDE = bIDE(:);
% bIDE, iIDE are indices into DE.
%
% interior undirected edges:
iiUE = find(DE(iIDE,1)<DE(iIDE,2));
iiRE = find(DE(iIDE,1)>DE(iIDE,2));
% iiUE is indices into iIDE
DE(iIDE(iiUE),:) % list of undirected edges, containing vertice indices.
DE(iIDE(iiRE),:) % list of undirected edges, containing vertice indices.
assert(size(iiUE,1)==size(iiRE,1));

if false
    % this is *wrong**
    [CC,~,iiUE_opposite] = intersect(DE(iIDE(iiUE),:),DE(iIDE(iiRE),[2,1]),'rows','stable');
    assert(size(iiUE,1)==size(CC,1));
    % such that DE(iIDE(iiUE),:)== DE(iIDE(iiRE),[2,1])(iiUE_opposite,:) as
    % asserted
    DD = DE(iIDE(iiRE),[2,1]);
    assert(norm(DE(iIDE(iiUE),:)-DD(iiUE_opposite,:))<1e-8);
else
    [CC,~,iiUE_opposite] = intersect(DE(iIDE(iiRE),[2,1]),DE(iIDE(iiUE),:),'rows','stable');
    assert(size(iiUE,1)==size(CC,1));
    % asserted
    DD = DE(iIDE(iiUE),:);
    assert(norm(DD(iiUE_opposite,:)-DE(iIDE(iiRE),[2,1]))<1e-8);
end 
    
%
iF_DE(iIDE(iiUE));
iF_DE(iIDE(iiRE));
iiF_DE(iIDE(iiUE));
iiF_DE(iIDE(iiRE));
% mid-point value constraints
% one per interior edge
t = (1:size(iiUE,1))';
% indices into DE
cIndex = iIDE(iiUE);
pIndex = mod(cIndex-f-1,3*f)+1;
nIndex = mod(cIndex+f-1,3*f)+1;
cIndex2 = iIDE(iiRE);
pIndex2 = mod(cIndex2-f-1,3*f)+1;
nIndex2 = mod(cIndex2+f-1,3*f)+1;
%
tb = (1:size(bIDE,1))';
% indices into DE
cBIndex = bIDE;
pBIndex = mod(cBIndex-f-1,3*f)+1;
nBIndex = mod(cBIndex+f-1,3*f)+1;
assert(norm(iF_DE(cBIndex)-iF_DE(pBIndex))+norm(iF_DE(cBIndex)-iF_DE(nBIndex))<1e-8)
assert(norm(iF_DE(cIndex)-iF_DE(pIndex))+norm(iF_DE(cIndex)-iF_DE(nIndex))<1e-8)
assert(norm(iF_DE(cIndex2)-iF_DE(pIndex2))+norm(iF_DE(cIndex2)-iF_DE(nIndex2))<1e-8)
assert(max(DE(cIndex,1)-DE(cIndex,2))<0)
assert(min(DE(cIndex2,1)-DE(cIndex2,2))>0)
%
sI = [t;t;t;...
    iiUE_opposite;iiUE_opposite;iiUE_opposite];
sJ = [cIndex;pIndex;nIndex;...
    cIndex2;pIndex2;nIndex2];
sV = [... *no* minus sign for the second triangle since normal derivatives are opposite at adjacent triangles. 
    -1./diag(EH(iF_DE(cIndex),iiF_DE(cIndex)));...
    diag(Ecos(iF_DE(nIndex),iiF_DE(nIndex)))./diag(EH(iF_DE(pIndex),iiF_DE(pIndex)));...
    diag(Ecos(iF_DE(pIndex),iiF_DE(pIndex)))./diag(EH(iF_DE(nIndex),iiF_DE(nIndex)));...
    -1./diag(EH(iF_DE(cIndex2),iiF_DE(cIndex2)));...
    diag(Ecos(iF_DE(nIndex2),iiF_DE(nIndex2)))./diag(EH(iF_DE(pIndex2),iiF_DE(pIndex2)));...
    diag(Ecos(iF_DE(pIndex2),iiF_DE(pIndex2)))./diag(EH(iF_DE(nIndex2),iiF_DE(nIndex2)));...
    ];
% sV = [EL(iF_DE(iIDE(iiUE));iiF_DE(iIDE(iiUE)))];

% each row for an undirected edge
% each col for an directed edge, total 3f of them.
GL = sparse(sI,sJ,sV,size(iiUE,1),3*f);
%
sbI = [tb;tb;tb;];
sbJ = [cBIndex;pBIndex;nBIndex;];
sbV = [... *no* minus sign for the second triangle. 
    2*(diag(Ecot(iF_DE(nBIndex),iiF_DE(nBIndex))).*diag(Ecot(iF_DE(pBIndex),iiF_DE(pBIndex))))./diag(EL(iF_DE(pBIndex),iiF_DE(pBIndex))).^2;...
    - bsxfun(@times, diag(Ecot(iF_DE(pBIndex),iiF_DE(pBIndex))),1./dblA(iF_DE(cBIndex),:));...
    - bsxfun(@times, diag(Ecot(iF_DE(nBIndex),iiF_DE(nBIndex))),1./dblA(iF_DE(cBIndex)));...
    ];
%
BQ = sparse(sbI,sbJ,sbV,size(bIDE,1),3*f);
%
sV = [... *no* minus sign for the second triangle.
    0.5./diag(EH(iF_DE(cIndex),iiF_DE(cIndex)));...
    -0.5./diag(EH(iF_DE(cIndex),iiF_DE(cIndex)));...
    -0.5./diag(EH(iF_DE(cIndex),iiF_DE(cIndex)));...
    0.5./diag(EH(iF_DE(cIndex2),iiF_DE(cIndex2)));...
    -0.5./diag(EH(iF_DE(cIndex2),iiF_DE(cIndex2)));...
    -0.5./diag(EH(iF_DE(cIndex2),iiF_DE(cIndex2)));...
    ];
%
GQ = sparse(sI,sJ,sV,size(iiUE,1),3*f);
%
sV = [... minus sign for the second triangle. 
    zeros(size(cIndex,1),1);...
    0.5*ones(size(cIndex,1),1);...
    0.5*ones(size(cIndex,1),1);...
    -zeros(size(cIndex2,1),1);...
    -0.5*ones(size(cIndex2,1),1);...
    -0.5*ones(size(cIndex2,1),1);...
    ];
%
VL = sparse(sI,sJ,sV,size(iiUE,1),3*f);
%
sV = [... minus sign for the second triangle. 
    0.25*ones(size(cIndex,1),1);...
    zeros(size(cIndex,1),1);...
    zeros(size(cIndex,1),1);...
    -0.25*ones(size(cIndex2,1),1);...
    zeros(size(cIndex2,1),1);...
    zeros(size(cIndex2,1),1);...
    ];
%
VQ = sparse(sI,sJ,sV,size(iiUE,1),3*f);
%
BL = sparse([],[],[],size(bIDE,1),3*f);
%
%
Ge = zeros(size(iiUE,1),1);% interior undirected edges
Be = zeros(size(bIDE,1),1);% boundary edges: second order boundary condition
%
L = [GL;VL;BL];
Q = [GQ;VQ;BQ];
%%
r = struct();
r.L = L;
r.Q = Q;

%%
if false
%%
h = max(max(V)-min(V));
u = V(:,1).*(V(:,1)/h).*3 + V(:,2).*(V(:,2)/h).^3 - V(:,1).*(V(:,1)/h).^5;
%u = V(:,2);
cl = [u(F(:,1),:),u(F(:,2),:),u(F(:,3),:)];
cq = Q\([Ge;Ge;Be]-L*cl(:));
cq = reshape(cq,[size(cq,1)/3,3]);
norm([Ge;Ge;Be]-L*cl(:)-Q*cq(:))
%%
residual = ([Ge;Ge;Be]-L*cl(:));
norm(residual)
end