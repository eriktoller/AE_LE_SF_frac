function [Omega] = Omega_dp(chi,b,beta,alpha,m)
% Set the initial value to zero
Omega = 0;
% Check is it can be calcuated with the far-field approximation
if chi*conj(chi) < 1+(1e-5)^2
    for jj = 0:m-1
        nn = jj+1;
        Omega = Omega + b(nn)*F_func(chi,beta,jj);
    end
else
    for jj = 1:m
        Omega = Omega + b(jj)*F_far(chi,alpha(jj,:));
    end
end
Omega = Omega/(2*pi);
end

