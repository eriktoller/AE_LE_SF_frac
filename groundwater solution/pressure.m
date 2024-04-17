function [ p ] = pressure( z,h,rho,g,Omega_func, phi_func )
% Get the value for Omega
Omega = Omega_func(z);
% Get the value for phi
phi = phi_func(real(Omega));
% Calclaute the pressure
h = h + imag(z);
p = (phi - h)*rho*g; 
end