function [ Q ] = solve_Q( z_ref,Phi0,Phi_tunnel,Phi_gw_level,AM,Ne,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,rt,zt,t,z1,z2,c )
M = length(zt);
M1 = length(z1);
Q0=zeros(1,M+M1);
KN=zeros(M+M1+1,1);
C = 0;
for jj=1:M
    temp=real(Cauchy_integral(Ne,0,0,@(Z)z_from_bigz(Z,rt(jj),zt(jj)),...
        @(z)Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q0,z1,z2,c,0,0,0,0)));
    KN(jj,1)=Phi_tunnel(jj)-temp(1);
end

for jj=1:M1
    temp=real(Cauchy_integral(Ne,0,0,@(chi)z_from_chi(chi,z1(jj),z2(jj)),...
        @(z)Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q0,z1,z2,c,0,0,0,0)));
    KN(M+jj,1)=Phi_gw_level(jj)-temp(1);
end

KN(M+M1+1,1)=Phi0-real(Omega_total(z_ref,Phi0*0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q0,z1,z2,c,0,0,0,0));

Q=AM\KN;


end

