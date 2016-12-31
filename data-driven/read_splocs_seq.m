function [VFs,F] = read_splocs_seq(fname)
%%
%fname = 'face.h5'

h5disp(fname);
verts = h5read(fname,'/verts');
tris = h5read(fname,'/tris');
% mean = h5read(fname,'/mean'); % no such field.
% scale = h5read(fname,'/scale'); % no such field.

VFs = permute(verts,[2,1,3]);
% the inverse mapping is: verts = permute(VFs,[2 1 3]);
F = permute(tris,[2 1]);
% the inverse mapping is: tris = permute(F,[2 1]);

assert(size(VFs,2)==3);
assert(size(F,2)==3);
% assert(size(mean,1)==1&&size(mean,2)==3);
% assert(length(scale)==1);

VFs = double(VFs);
F = double(F) + 1; % change from 0-index to 1-index.