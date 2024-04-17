function [a] = solve_acoef(theta,W_uni,z1b,z2b,b,z1a,z2a,a_in,beta,alpha,m,N,A,pos,rt,zt,t,Q,z1gw,z2gw,c,m_not)
% Get the evaluation points
mu = angle(z2a(m_not)-z1a(m_not));
chi = exp(1i*theta);
z = z_from_chi(chi,z1a(m_not),z2a(m_not));

% Calcualte the value for W at the evalutaion points
w = zeros(N,1);
for ii = 1:N
    w(ii) = imag(W_total(z(ii),W_uni,z1b,z2b,b,z1a,z2a,a_in,beta,alpha,m,rt,zt,t,Q,z1gw,z2gw,c,0,m_not,0,0)*exp(1i*mu));
end


% Solve for the coefficients a
a = A\w;
end

