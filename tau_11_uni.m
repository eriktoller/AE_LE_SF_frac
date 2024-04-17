function [tau_11] = tau_11_uni(z,W_uni,gamma_w,k,nu,kappa)
%Tau 11 for uniform seepage flow in a half-space
%   Detailed explanation goes here
tau_11 = 2i*W_uni*gamma_w/k*(1-2*nu)/(kappa+1)*(z-conj(z));
end

