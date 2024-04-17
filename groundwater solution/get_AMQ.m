function [ AM ] = get_AMQ(N,rt,zt,z1,z2,z_ref )
M = length(zt);
M1 = length(z1);
for ii = 1:M
    Z_far(ii) = bigz_from_z(z_ref,rt(ii),zt(ii));
end
for ii = 1:M1
    chi_far(ii) = chi_from_z(z_ref,z1(ii),z2(ii));
end
AM=zeros(M+M1+1,M+M1+1);
for jj=1:M
    for mm=1:M
        a=real(Cauchy_integral(N,0,0,@(Z)z_from_bigz(Z,rt(jj),zt(jj)),...
            @(z)log(bigz_from_z(z,rt(mm),zt(mm))/abs(Z_far(mm)))/(2*pi)));
        AM(jj,mm)=a(1);
    end
    for mm=1:M1
        a=real(Cauchy_integral(N,0,0,@(Z)z_from_bigz(Z,rt(jj),zt(jj)),...
            @(z)log(chi_from_z(z,z1(mm),z2(mm))/abs(chi_far(mm)))/(2*pi)));
        AM(jj,M+mm)=a(1);
    end
    AM(jj,M+M1+1)=1;
end
for jj=1:M1
    for mm=1:M
        a=real(Cauchy_integral(N,0,0,@(chi)z_from_chi(chi,z1(jj),z2(jj)),...
            @(z)log(bigz_from_z(z,rt(mm),zt(mm))/abs(Z_far(mm)))/(2*pi)));
        AM(M+jj,mm)=a(1);
    end
    for mm=1:M1
        a=real(Cauchy_integral(N,0,0,@(chi)z_from_chi(chi,z1(jj),z2(jj)),...
            @(z)log(chi_from_z(z,z1(mm),z2(mm))/abs(chi_far(mm)))/(2*pi)));
        AM(M+jj,M+mm)=a(1);
    end
    AM(M+jj,M+M1+1)=1;
end

for jj=1:1
    for mm=1:M
            a=real(log(bigz_from_z(z_ref,rt(mm),zt(mm))/abs(Z_far(mm)))/(2*pi));
        AM(M+M1+jj,mm)=a(1);
    end
    for mm=1:M1
            a=real(log(chi_from_z(z_ref,z1(mm),z2(mm))/abs(chi_far(mm)))/(2*pi));
        AM(M+M1+jj,mm)=a(1);
    end
    AM(M+M1+jj,M+M1+1)=1;
end
end

