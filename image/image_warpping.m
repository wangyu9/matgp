function [img_out] = image_warpping(img_in,V,F,V2)

V2=V2*1.0;%image factor
%V2 = bsxfun(@minus,V2,min(V2));
new_size = max(ceil(V2));
mm = new_size(1);
nn = new_size(2);
[aa,bb,~]= size(img_in);

V(:,1) = aa-V(:,1);
V2(:,1) = mm-V2(:,1);

% V(:,2) = bb+1-V(:,2);
% V2(:,2) = nn+1-V2(:,2);

img_out = 255*ones(new_size(1),new_size(2),3,'uint8'); 
%imshow(img_in);
%img_out = img_in;
for i=1:1:size(F,1)
    V2(F(i,:)',:);
    mmin = min(V2(F(i,:)',:));
    mmax = max(V2(F(i,:)',:));
    rect = span_2D(mmin,mmax);
    %rect(:,1) = aa+1-rect(:,1);
    %rect(:,2) = bb - rect(:,2);
    P = [V2(F(i,:)',:)';[1,1,1];];
    lamada = P\[rect';ones(1,size(rect,1));];
    lamada = lamada';
    inside_list = find(lamada(:,1)>=0&lamada(:,2)>=0&lamada(:,3)>=0);
    %inside_list = (1:size(lamada,1))';
    %outside_list = find(~sparse(1,zero_index,true,1,n));
    xx = lamada(inside_list,:)*[V(F(i,1),1);V(F(i,2),1);V(F(i,3),1)];
    yy = lamada(inside_list,:)*[V(F(i,1),2);V(F(i,2),2);V(F(i,3),2)];
    XX = round(xx);
    YY = round(yy);
    [XX,YY] = tune(XX,YY,aa,bb);
    XX_dst = rect(inside_list,1);
    YY_dst = rect(inside_list,2);
    [XX_dst,YY_dst] = tune(XX_dst,YY_dst,mm,nn);
    img_out(XX_dst,YY_dst,:) = img_in(XX,YY,:);
%     img_out(XX,YY,1) = img_in(XX,YY,1);
%     img_out(XX,YY,2) = img_in(XX,YY,2);
%     img_out(XX,YY,3) = img_in(XX,YY,3);
end
imshow(img_out);
hold on;
for i=1:size(F,1)
    line(V2(F(i,1:2),2),V2(F(i,1:2),1));
    line(V2(F(i,2:3),2),V2(F(i,2:3),1));
    line(V2(F(i,[1,3]),2),V2(F(i,[1,3]),1));
end

end

function [XX,YY] = tune(XX,YY,aa,bb)
    XX(find(XX<=0))=1;
    YY(find(YY<=0))=1;
    XX(find(XX>aa))=aa;
    YY(find(YY>bb))=bb;
end

function [list] = span_2D(mmin,mmax)
    mmin = floor(mmin);
    mmax = ceil(mmax);
    disp = mmax-mmin;
    m = disp(1);
    n = disp(2);
%     width = (max(disp(:,1))+1);
%     list = bsxfun(@plus,0:1:disp(:,1),(0:1:disp(:,2))'*width);
%     list = list+mmin(:,1)+width;
    %list(:,1) = 
    X = ones(n+1,1)*(mmin(1):mmax(1));
    Y = (ones(m+1,1)*(mmin(2):mmax(2)))';
    list = [X(:),Y(:)];
end