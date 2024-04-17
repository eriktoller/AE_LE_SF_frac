function [W] = W_dp(chi,z1,z2,b,beta,alpha,m)
% Calcuate the prerequisits
W = 0;
mu = angle(z2-z1);
L = sqrt((z2-z1)*conj(z2-z1));
Theta = 4*chi^2/( L*(chi^2-1) )*exp(-1i*mu);
% Check is the far field expansion should be used
if chi*conj(chi) < 1+(1e-5)^2
    for jj = 0:m-1
        nn = jj+1;
        W = W - b(nn)*F_prim_func(chi,beta,jj);
    end
else
    for jj = 1:m
        W = W - b(jj)*F_prim_far(chi,alpha(jj,:));
    end
end
% Add constants to W
W = W*Theta/(2*pi);
end

