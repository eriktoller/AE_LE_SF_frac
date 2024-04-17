function [b] = solve_bcoef(theta,W_uni,z1b,z2b,b_in,z1a,z2a,a,beta,alpha,m,N,A,pos,rt,zt,t,Q,z1gw,z2gw,c,m_not)
% Get the evaluation points
mu = angle(z2b(m_not)-z1b(m_not));
chi = exp(1i*theta);
z = z_from_chi(chi,z1b(m_not),z2b(m_not));

% Calcualte the value for W at the evalutaion points
w = zeros(N,1);
for ii = 1:N
    w(ii) = real(W_total(z(ii),W_uni,z1b,z2b,b_in,z1a,z2a,a,beta,alpha,m,rt,zt,t,Q,z1gw,z2gw,c,m_not,0,0,0)*exp(1i*mu));
end
n = cos((0:m-1)*pi);
nab = repmat(n,size(b_in,1),1);
bb = b_in;
bb(m_not,:) = [];
error_ab1 = sum(bb(pos(1,:),:).*nab(pos(1,:),:),'all') - sum(bb(pos(2,:),:),'all');
error_ab2 = sum(bb(pos(3,:),:).*nab(pos(3,:),:),'all') - sum(bb(pos(4,:),:),'all');
w(end+1) = -error_ab1;
w(end+1) = -error_ab2;

%[ w1 ] = w_bulider( m,w(1:end-2),m_not,z1b,z2b,b_in );

% Solve for the coefficients b
b = A\w;
end

