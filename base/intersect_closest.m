function [d,I] = intersect_closest(A,B)

A = A';
B = B';

Dist = zeros(2,size(A,2));
for p = 1:size(A,2)
 [m,q] = min(sum(abs(bsxfun(@minus,B,A(:,p))),1));
 Dist(:,p) = [m;q];
end

d = Dist(1,:)';
I = Dist(2,:)';