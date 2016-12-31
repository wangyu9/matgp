function [t] = render_animation(VFs,F)

if(size(VFs,3)==1)
    assert(mod(size(VFs,1),3)==0);
    VFs = reshape(VFs,[size(VFs,1)/3,3,size(VFs,2)]);
end

assert(size(VFs,2)==3);

nf = size(VFs,3);

hold off;
axis manual;

[t,~,~] = render_mesh(VFs(:,:,1),F,'view',[90 90]);
%%
for i=1:nf
   set(t,'Vertices',VFs(:,:,i));
   drawnow
end