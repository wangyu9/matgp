function [] = writeH5MAT(fname,W)

verbose = false;

h5create(fname,'/H5MATRIXDENSE',size(W));

h5write(fname,'/H5MATRIXDENSE',W,[1 1],size(W));

if(verbose)
    h5disp(fname);
end