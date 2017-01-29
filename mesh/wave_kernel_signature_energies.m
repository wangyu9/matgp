function [log_Es] = wave_kernel_signature_energies(DD,N)
% N is the number of samples.
log_E = log(max(abs(diag(DD)),1e-6))';
log_Es = linspace(log_E(2),(max(log_E))/1.02,N);