function [Ts] = heat_kernel_signature_times(DD,N)
% N is the number of samples.

t_min = 4*log(10)/DD(end,end);
t_max = 4*log(10)/DD(2,2);

assert(t_min<t_max);

Ts = exp(linspace(log(t_min),log(t_max),N));