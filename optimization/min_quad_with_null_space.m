function [Z] = min_quad_with_null_space(A,B,zero_index,Aeq,Beq)

n = size(A,1);

Z = zeros(n,1);

none_zero_index = find(~sparse(1,zero_index,true,1,n));

Ann = A(none_zero_index,none_zero_index);

Bn = B(none_zero_index);

% some constraints should be removed since they becomes all zero
Aeq_temp = Aeq(:,none_zero_index);
non_zero_rows = find(sum(abs(Aeq_temp),2)~=0);
min(sum(abs(Aeq_temp),2))

Aeq_new = Aeq_temp(non_zero_rows,:);
Beq_new = Beq(non_zero_rows);

sum(sum(Ann~=0));

%[W] = min_quad_with_fixed(Ann,Bn,[],[],Aeq_new,Beq);

[SLeft SRight] = spspaces(Aeq_new,2,1e-22);
Stprime = SRight{1};
St = Stprime(:,SRight{3});
max(max(abs(Aeq_new*St))) % double check the max error of null space decomposition

W0 = mldivide(Aeq_new,Beq_new);

[Y] = min_quad_with_fixed(St'*Ann*St,St'*Bn+2*St'*Ann*W0,[],[],[],[]);

W = W0+St*Y;

%         for round=1:2
%             slots = 4;
%             for i=0:(slots-1);
%                 known = cols;
%                 Y = zeros(size(known,1),1);
%                 known_id = ((i*total/slots+1):(i+1)*total/slots)';
%                 known = [known;known_id];
%                 known = unique(known);
%                 Y = W(known);
%                 [W] = min_quad_with_fixed_layer(Q,zeros(n*m,1),known,Y,Aeq,Beq);
%                 %[W,F,Lambda,Lambda_known] = min_quad_with_fixed(Q,zeros(n*m,1),known,Y,Aeq,Beq);
%             end
%         end

Z(none_zero_index) = W;