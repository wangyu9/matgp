function [HKS] = heat_kernel_signature_dirac(Vs,DD,Ts,varargin)
% VV is of size n*k
% EV is of size n*4*k
% DD is of size k*k

    n = size(EV,1);
    dim = size(EV,2);
    assert(dim==3);
 
% there is no need to design new ones for the diracian case. 

error(['call heat_kernel_signature directly!']);

HKS = heat_kernel_signature_dirac(Vs,DD,Ts,varargin);