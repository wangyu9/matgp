function [t1] = compare_animation(VF1s,F1,VF2s,F2)

if(size(VF1s,3)==1)
    assert(mod(size(VF1s,1),3)==0);
    VF1s = reshape(VF1s,[size(VF1s,1)/3,3,size(VF1s,2)]);
end

if(size(VF2s,3)==1)
    assert(mod(size(VF2s,1),3)==0);
    VF2s = reshape(VF2s,[size(VF2s,1)/3,3,size(VF2s,2)]);
end

assert(size(VF1s,2)==3);
assert(size(VF2s,2)==3);
assert(size(VF1s,3)==size(VF2s,3));

nf = size(VF1s,3);

hold off;
axis manual;


subplot(1,3,1);
[t1,~,~] = render_mesh(VF1s(:,:,1),F1,'view',[90 90]);

subplot(1,3,2);
[t2,~,~] = render_mesh(VF2s(:,:,1),F2,'view',[90 90]);

subplot(1,3,3);
[t3,~,~] = render_mesh(VF1s(:,:,1),F1,'view',[90 90]);
colormap(my_colormap('weights-neg'));
%%
for i=1:nf
   set(t1,'Vertices',VF1s(:,:,i));
   set(t2,'Vertices',VF2s(:,:,i));
   set(t3,'Vertices',VF1s(:,:,i));
   
   % set t3's color
    dif = sum( (VF1s(:,:,i) - VF2s(:,:,i)).^2, 2);
    dif = dif/max(abs(dif));
    
%    caxis([-0.2,1]);
%    colorbar
    set(t3,'CData',dif);
   
   drawnow
end