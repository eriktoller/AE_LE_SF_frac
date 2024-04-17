function [w] = w_disp(z,z1,z2,L,m,beta,mu,sigma_11,z1b,z2b,b,alpha,W_uni,gamma_w,k,kappa,G,nu)
wf = complex(0,0);
for ii = 1:length(z1)
    wf = wf + w_frac(z,z1(ii),z2(ii),L(ii),m,beta(ii,:),mu(ii),kappa,G);
end
wp = w_uni(z,sigma_11,kappa,G)*0;

B2 = z*z*W_uni;
for ii = 1:length(z1b)
    chi = chi_from_z(z,z1b(ii),z2b(ii));
    for kk = 1:length(b)
        for jj = 1:length(alpha(kk,:))
            F = alpha(kk,jj)*chi^(-(jj-1))/(1-jj) + chi^(-(jj+1))/(jj+1);
        end
        B2 = B2 + b(kk)*F;
    end
end
B2 = B2 * gamma_w/(2*k);

w = 1/(4*G)*(wf + wp + 2*(1-2*nu)*conj(B2));

end

