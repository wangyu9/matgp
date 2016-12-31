function [t] = draw_weights(V,F,W,varargin)

assert(size(V,1)==size(W,1));
assert(size(W,2)==1);

% input parser.
C = [];

nvar = length(varargin);
ii=1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'C'))
       C = varargin{ii+1};
       ii = ii + 1;
   end
   ii = ii + 1;
end


%%
assert(max(-Inf,0)==0);
%cf = @(x) ( (x>=0).*max(log(abs(x)),-5)+(x<0).*(-10-max(log(abs(x)),-5)) + 5)./5;
%W = cf(V(:,3));

%W = log(abs(V(:,3)));
%W(W<-3)=0;

%%


t = trisurf(F,V(:,1),V(:,2),0*V(:,1),'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
hold on;

if(size(C,1)>0)
    % for me each coordinate corresponds to a point in C
    scatter(C(:,1),C(:,2),'ok','MarkerFaceColor','y','LineWidth',3,'SizeData',100);
end

hold off;
axis equal;
view(2);


exjet = my_colormap('weights-neg'); % remerber to have caxis([-0.2,1]);
colormap(exjet);%colormap('jet');
if(false)
    colorbar
end   
    %title('Vector Field  ','FontSize',20);


set(t,'CData',W);
drawnow;

caxis([-0.2,1]);
grid off;
axis off;
