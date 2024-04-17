function [tau_12] = tau_12_uni(z,W_uni,gamma_w,k,kappa)
%Tau 12 for uniform seepage flow in a half-space
%   Detailed explanation goes here
tau_12 = -2i*W_uni*gamma_w/k*1/(kappa+1)*(z-conj(z));
end

