function [tau_12] = tau_12_uni(z,W_uni,gamma_w,k,kappa)
tau_12 = -2i*W_uni*gamma_w/k*1/(kappa+1)*(z-conj(z));
end

