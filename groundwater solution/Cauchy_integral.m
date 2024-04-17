function [ a ] = Cauchy_integral( N_in,m,m_max,z_of_chi,Omega_of_z )
if N_in<2*m
    N=2*m;
else
    N=N_in;
end
deltheta=2*pi/N;
theta_0=0.5*deltheta;
Int=zeros(N,m_max+1);
a=zeros(1,m_max+1);
for nu=1:N
    n=nu-1;
    theta=theta_0+n*deltheta;
    chi=exp(1i*theta);
    z=z_of_chi(chi);
    Omega=Omega_of_z(z);
    for j=1:m+1
        mu=j-1;
        Int(nu,j)=real(Omega)*exp(-1i*mu*theta);
    end
end
for j=1:m+1
   a(j)=0;
   for n=1:N
     a(j)=a(j)+Int(n,j);
   end
   a(j)=2*a(j)/N;
end
a(1)=.5*a(1);
end

