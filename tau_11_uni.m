function [tau_11] = tau_11_uni(z,W_uni,gamma_w,k,nu,kappa)

tau_11 = 2i*W_uni*gamma_w/k*(1-2*nu)/(kappa+1)*(z-conj(z));
end

