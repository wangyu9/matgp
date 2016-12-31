function [A_data,B_data] = QuadMeshSequence(VFs,CFs)

% the output is the quadratic objective such that Error = w'Aw + w'B + C
% C is to be implemented.

n = size(VFs,1);
d = size(VFs,2);
nf = size(VFs,3);

m = size(CFs,1);
assert(d==size(CFs,2));
assert(nf==size(CFs,3));

%A_data = spalloc(n*m,m*m*n);
%B_data = spalloc(n*m,1,m*n);

IA = zeros(n,m*m);
JA = zeros(n,m*m);
KA = zeros(n,m*m);

IB = zeros(n,m);
JB = zeros(n,m);
KB = zeros(n,m);


for j=1:1:n % for all vertices
   
   sumKA = zeros(m*m,1);
   sumKB = zeros(m,1); 
    
   for k =1:1:d % for all dimensions

       % the following codes is doing this in nature:
       % A_data = A_data + c_k * c_k';
       % B_data = B_data - 2 * v_ijk * c_k;

       for i=1:1:nf % for all frames
           v_ijk = VFs(j,k,i);
           sumKA = sumKA + reshape( CFs(:,d,i) * CFs(:,d,i)', [m*m,1]);
           sumKB = sumKB - 2 * v_ijk * CFs(:,d,i);
       end
   end
   
   IA(j,:) = reshape( repmat((0:m-1)'*n+j,1,m), [m*m,1]);
   JA(j,:) = reshape( repmat((0:m-1)'*n+j,1,m)', [m*m,1]);
   KA(j,:) = sumKA;

   IB(j,:) = (0:m-1)'*n+j;
   JB(j,:) = ones(m,1);
   KB(j,:) = sumKB;
       
end
    
B_data = sparse(IB(:),JB(:),KB(:),n*m,1);
A_data = sparse(IA(:),JA(:),KA(:),n*m,n*m);