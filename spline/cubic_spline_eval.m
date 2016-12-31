function [V,I,M] = cubic_spline_eval(tC,tV,CM,C)

m = size(CM,1)+1;
assert(m==size(C,1));
assert(size(tV,2)==1);
n = size(tV,1);
V = zeros(n,size(C,2));
I = zeros(n,1);
% do not need this assert(max(tV)<=max(tC));
for i=1:1:n
    left = 1;
    for j=1:1:m-1
        if(tV(i)<=tC(j+1))
            left = j;
            break;
        else
            left = j+1;
        end
    end
    I(i) = left;
    
    V(i,:) = [1,tV(i),tV(i)^2,tV(i)^3] * CM{left} * C;
end


M = zeros(n,m);
for i=1:1:m-1
    J = find(I==i);
    M(J,:) = [tV(J,:).^0, tV(J,:).^1, tV(J,:).^2, tV(J,:).^3] * CM{i};
end