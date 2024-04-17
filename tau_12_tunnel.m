function [tau_12] = tau_12_tunnel(z,rt,zt,a,b,g,rho,rho_w,k,Omega_w_func)
%Tau 11 for uniform seepage flow in a half-space
%   Detailed explanation goes here
Z = bigz_from_z(z,rt,zt);
if Z*conj(Z) < 1-1e-12
    tau_12 = complex(NaN,NaN);
else
dPhi = 0;
Om = Omega(z,g,rho,rho_w,k,Omega_w_func);
for ii = 1:length(a)
    nn = ii-1;
    dPhi = dPhi + ( b(ii) - a(ii))*Z^(-nn);
end
tau_12 = -dPhi - conj(dPhi) + 0*(conj(Om) + Om);


end

end

