function [Omega] = Omega_db(chi,a,beta,alpha,m)
% Set the initial value to zero
Omega = 0;
% Check if Omega can be calculated with the far-field approximation
if chi*conj(chi) < 1+(1e-5)^2
    for jj = 0:m-1
        nn = jj+1;
        Omega = Omega + a(nn)*F_func(chi,beta,jj);
    end
else
    for jj = 1:m
        Omega = Omega + a(jj)*F_far(chi,alpha(jj,:));
    end
end
Omega = Omega/(2i*pi);
end

