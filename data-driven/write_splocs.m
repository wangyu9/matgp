function [] = write_splocs(fname,Ws,V,F)

verbose = false;

if(size(Ws,3)==1)
    assert(mod(size(Ws,1),3)==0);
    Ws = reshape(Ws,[size(Ws,1)/3,3,size(Ws,2)]);
end

assert(size(Ws,1)==size(V,1));
assert(size(Ws,2)==3);
assert(size(V,2)==3);
assert(size(F,2)==3);

decom = permute(Ws,[2 1 3]);
tris = permute(F,[2 1]);
verts = permute(V,[2 1]);

decom = bsxfun(@plus,decom,verts);

tris = tris - 1; % made it from 1-indexed to 0-indexed.

for i=1:1:size(decom,3)
   h5create(fname, ['/comp',sprintf('%03d',i-1)], [size(decom,1),size(decom,2)]);
end

for i=1:1:size(decom,3)
   h5write(fname, ['/comp',sprintf('%03d',i-1)], permute(decom(:,:,i),[1 2]), [1 1], [size(decom,1),size(decom,2)]);
end

h5create(fname,'/default',size(verts));
h5create(fname,'/tris',size(tris),'Datatype','int32'); 

h5write(fname,'/default',verts,[1 1],size(verts));
h5write(fname,'/tris',int32(tris),[1 1],size(tris));

%h5write(fname,'/verts',single(decom),[1 1 1],size(decom));
%h5write(fname,'/tris',int32(tris),[1 1],size(tris));
if(verbose)
    h5disp(fname);
end