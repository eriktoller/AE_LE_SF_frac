function [tau_11] = tau_11_tunnel(z,rt,zt,a,b)

Z = bigz_from_z(z,rt,zt);
if Z*conj(Z) < 1-1e-12
    tau_11 = complex(NaN,NaN);
else
ddPhi = 0;
dpsi = 0;
for ii = 1:length(a)
    nn = ii-1;
    ddPhi = ddPhi + nn*(b(ii) - a(ii))*Z^(-(nn+1));
    dpsi = dpsi + a(ii)*Z^(-(nn+2));
end
ddPhi = -(1/rt) * ddPhi;
dpsi = 2*dpsi;
tau_11 = rt*(conj(Z)-1/Z)*ddPhi - dpsi;
end

end

