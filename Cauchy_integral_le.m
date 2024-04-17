function [ aa ] = Cauchy_integral_le( ReIm,N_in,m,m_max,z_of_chi,g,rho,rho_w,k,Omega_w_func,T_func,m_not_t )
if N_in<2*m
    N=2*m;
else
    N=N_in;
end
deltheta=2*pi/N;
theta_0=0.5*deltheta;
Int=zeros(N,m_max+1);
aa=zeros(1,m_max+1);
for nu=1:N
    n=nu-1;
    theta=theta_0+n*deltheta;
    ang = theta + pi/2;
    chi=exp(1i*theta);
    z=z_of_chi(chi);
    T = T_func(z,ang,0,m_not_t);
    if ReIm
        for j=1:m+1
            mu=j-1;
            Om = Omega(z,g,rho,rho_w,k,Omega_w_func);
            Int(nu,j)=real(T)*exp(-1i*mu*theta);
        end
    else
        for j=1:m+1
            mu=j-1;
            Om = Omega(z,g,rho,rho_w,k,Omega_w_func);
            Int(nu,j)=(-imag(T) - .5*(Om + conj(Om)))*exp(-1i*mu*theta);
        end
    end
end
for j=1:m+1
   aa(j)=0;
   for n=1:N
     aa(j)=aa(j)+Int(n,j);
   end
   aa(j)=2*aa(j)/N;
end
aa(1)=.5*aa(1);
end

