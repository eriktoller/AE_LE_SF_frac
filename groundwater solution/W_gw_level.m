function [W] = W_gw_level(z,z1,z2,a,Q)
chi = chi_from_z(z,z1,z2);
if chi*conj(chi) < (1-1e-6)
    W = complex(NaN,NaN);
else
    W = 0;
    for nn = 1:length(a)-1
        W = W + nn*a(nn+1)/chi^(nn+1);
    end
    
    % Add the discharge
    W = W - Q/(2*pi)*(1/chi);
end
mu = angle(z2 - z1);
L = real(sqrt((z2 - z1) * conj(z2 - z1)));
Theta = 4.0 * chi * chi / (L * (chi * chi - 1.0)) * exp(-1i * mu);
W = W*Theta;
end

