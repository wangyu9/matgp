function [] = write_splocs_seq(fname,VFs,F)

%,mean,scale
assert(size(VFs,2)==3);
assert(size(F,2)==3);
%assert(size(mean,1)==1&&size(mean,2)==3);
%assert(length(scale)==1);

%verts = zeros([size(VFs,2) size(VFs,1) size(VFs,3)]);
%for i=1:1:size(VFs,3)
%   verts(:,:,i) = VFs(:,:,3)'; 
%end
verts = permute(VFs,[2 1 3]);

%tris = F';%zeros([size(F,2) size(F,1)]);
tris = permute(F,[2 1]);
tris = tris - 1; % made it from 1-indexed to 0-indexed.

h5create(fname,'/verts',size(verts),'Datatype','single');%,'Deflate',4);
h5create(fname,'/tris',size(tris),'Datatype','int32'); 
%h5disp(fname);

%h5writeatt(fname,'/','mean',mean); %e.g. [1 2 3]
%h5writeatt(fname,'/','scale',scale); %e.g. 4.4

h5write(fname,'/verts',single(verts),[1 1 1],size(verts));
h5write(fname,'/tris',int32(tris),[1 1],size(tris));

h5disp(fname);