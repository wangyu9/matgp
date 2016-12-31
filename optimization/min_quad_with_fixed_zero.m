function [Z] = min_quad_with_fixed_zero(A,B,zero_index,Aeq,Beq)

n = size(A,1);

Z = zeros(n,1);

none_zero_index = find(~sparse(1,zero_index,true,1,n));

Ann = A(none_zero_index,none_zero_index);

Bn = B(none_zero_index);

Aeq_new = Aeq(:,none_zero_index);

sum(sum(Ann~=0));

[W] = min_quad_with_fixed(Ann,Bn,[],[],Aeq_new,Beq);

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