function [Vm,Tm,Fm] = merge_mesh(V1,T1,F1,V2,T2,F2)
% merge two meshes into a single mesh
% if need to display anymesh before merge, could do this:
% bsxfun(@plus,V2,[400,0,0])


Vm = [V1;V2;];

%s = max(max(T1));
s = size(V1,1);
if(length(s)~=1)
    s = 0;
end
Tm = [T1;T2+s];

%s = max(max(T1));
if(length(s)~=1)
    s = 0;
end
Fm = [F1;F2+s];

end