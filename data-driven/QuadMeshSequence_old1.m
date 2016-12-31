function [A_data,B_data] = QuadMeshSequence(mesh_frames,C)

% the output is the quadratic objective such that Error = w'Aw + w'B

n = size(mesh_frames{1},1);
d = size(mesh_frames{1},2);

nf = size(mesh_frames,1);
m = size(C{1},1);

assert(nf==size(C,1));
assert(d==size(C{1},2));

%A_data = spalloc(n*m,m*m*n);
%B_data = spalloc(n*m,1,m*n);

IA = zeros(n,d,m*m);
JA = zeros(n,d,m*m);
KA = zeros(n,d,m*m);

IB = zeros(n,d,m);
JB = zeros(n,d,m);
KB = zeros(n,d,m);


    for j=1:1:n % for all vertices
       for k =1:1:d % for all dimensions
           v_ijk = mesh_frames{i}(j,k);
           
           c_k = sparse((0:m-1)'*n+j,ones(m,1),C{i}(:,d),n*m,1);
           % c_k((0:m-1)'*n+j,:) = C{i}(:,d); % this is too slow
           % A_data = A_data + c_k * c_k';
           % B_data = B_data - 2 * v_ijk * c_k;
           
           sumKA = zeros(m*m,1);
           sumKB = zeros(m,1);
           
           for i=1:1:1%nf % for all frames
               sumKA = sumKA + reshape( C{i}(:,d) * C{i}(:,d)', [m*m,1]);
               sumKB = sumKB - 2 * v_ijk * C{i}(:,d);
           end
           
           IA(j,k,:) = reshape( repmat((0:m-1)'*n+j,1,m), [m*m,1]);
           JA(j,k,:) = reshape( repmat((0:m-1)'*n+j,1,m)', [m*m,1]);
           KA(j,k,:) = sumKA;
           
           IB(j,k,:) = (0:m-1)'*n+j;
           JB(j,k,:) = ones(m,1);
           KB(j,k,:) = sumKB;
       end
    end
B_data = sparse(IB(:),JB(:),KB(:),n*m,1);
A_data = sparse(IA(:),JA(:),KA(:),n*m,n*m);