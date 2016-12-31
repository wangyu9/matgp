%function [] = skinning_decomposition(VFs,F)

% This implements the paper: Fast and Efficient Skinning of Animated
% Meshes
%%
threshold = 10;
K = 50;
%%
dim = 3;
A = reshape( permute(VFs,[2,3,1]), [nf*dim,n] );


%%
B = [];
C = [];
R = A;
i = 0;
while(1)
    if(i>=K)
        break;
    end
    [~,max_norm_col] = max(sum(R.^2,1));% find the column with max squared norm.
    mi = R(:,max_norm_col);
    if(i>0)
        % stabilization
        mi = mi - B*(B'*mi);
    end
    % normalization
    mi = mi/norm(mi);
    % deflation
    R = R - mi*mi'*R;
    B = [B,mi];
    C = [C;mi'*R];
    i = i + 1;
end

%%
[uu,ss,vv] = svd(A,K);
image(log(abs( A-uu*ss*vv' )),'CDataMapping','scaled');
%%
image(log(abs( A-B*C )),'CDataMapping','scaled');
%%
image(log(abs( B'*B-eye(size(B,2)) )),'CDataMapping','scaled');
hold on;
color map;
colorbar;