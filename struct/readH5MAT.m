function [W] = readH5MAT(fname)

verbose = false;

W = h5read(fname,'/H5MATRIXDENSE');

if(verbose)
    h5disp(fname);
end