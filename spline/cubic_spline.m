function [CM] = cubic_spline(tC)



assert(size(tC,2)==1);
m = size(tC,1);

A = sparse(m,m);
B = sparse(m,m);
A = zeros(m,m);
B = zeros(m,m);
A(1,1) = 1;
A(m,m) = 1;

b = (2:m-1)';

A = A + sparse(b, b-1, 1/6*(tC(b)-tC(b-1)) ,m, m);
A = A + sparse(b, b, 1/3*(tC(b+1)-tC(b-1)), m, m);
A = A + sparse(b, b+1, 1/6*(tC(b+1)-tC(b)) ,m, m);
B = B + sparse(b, b-1, 1./(tC(b)-tC(b-1)) ,m, m);
B = B + sparse(b, b, -1./(tC(b)-tC(b-1))-1./(tC(b+1)-tC(b)) ,m, m);
B = B + sparse(b, b+1, 1./(tC(b+1)-tC(b)) ,m, m);

D2= A\B;

coeffs = zeros(4,m-1);
CM = cell(m-1,1);

for i=1:1:m-1
   X = zeros(4,4);
   X(1,:) = [1,tC(i),tC(i)^2,tC(i)^3];
   X(2,:) = [1,tC(i+1),tC(i+1)^2,tC(i+1)^3];
   X(3,:) = [0,0,2,6*tC(i)];
   X(4,:) = [0,0,2,6*tC(i+1)];
   R = inv(X);
   %coeffs(:,i); 
   CM{i}= zeros(4,m);
   CM{i} = CM{i} + R(:,3) * D2(i,:) + R(:,4) * D2(i+1,:);
   CM{i}(:,i) = CM{i}(:,i) + R(:,1);
   CM{i}(:,i+1) = CM{i}(:,i+1) + R(:,2);
end