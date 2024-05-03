function [Om] = Omega(z,g,rho,rho_w,k,Omega_w_func)
%Omega for body force
Omega_w = Omega_w_func(z);
Om = rho*g*(1i*z - rho_w/(rho*k)*Omega_w);
end

