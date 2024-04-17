function [beta,a,b,theta,p,z] = solv_coef(cond,z1,z2,L,m,N,mu,rt,zt,sigma_11,g,rho,rho_w,k,h,Omega_w_func,Om_not)
% Creating empty Matrices
T = zeros(N,1);
beta = zeros(length(z1),m);
a = zeros(length(zt),m+1);
b = zeros(length(zt),m+1);

% TEMP
p = ones(N,length(z1));
pg = zeros(N,length(z1));
Om = zeros(N,length(z1));

error = 1;
beta_old = beta;
a_old = a;
b_old = b;
NIT = 0;


% Steps for the path along the unit circle in the chi-plane
theta_0 = pi/N;
deltheta = (pi - 2*theta_0)/(N-1);

% Getting the A matrix and theta
A = zeros(N,m);
z = zeros(length(z1),m);
theta = zeros(N,1);
term = zeros(N,length(z1));
for ii = 1:N
    theta(ii) = theta_0 + ii*deltheta;
    for mm = 1:m
        A(ii,mm) = mm*sin( mm*theta(ii) );  
    end
    for jj = 1:length(z1)
        chi = exp(1i*theta(ii));
        z(jj,ii) = z_from_chi(chi,z1(jj),z2(jj));
        term(ii,jj) = -0.5*L(jj)*sin(theta(ii));
        p(ii,jj) = (pg(ii,jj) - pressure( z(jj,ii),h,rho_w,g,Omega_w_func, @(Phi) fi_from_Phi( Phi,k ) ) );
        if Om_not == 1
            Om(ii,jj) = Omega(z(jj,ii),g,rho,0*rho_w,k,Omega_w_func);
        else
            Om(ii,jj) = Omega(z(jj,ii),g,rho,rho_w,k,Omega_w_func);
        end
    end
end

cnt = 0;
error_old = error;
str = '     Error    Iteration';
disp(str)

% iteration solver
while error > cond && NIT < 300
    % Caculating the T from all other elements
    for ii = 1:length(z1)
        
        for jj = 1:N
            T(jj) = T_total(z(ii,jj),z1,z2,L,m,beta,mu,rt,zt,a,b,mu(ii),sigma_11,g,rho,rho_w,k,Omega_w_func ,ii,0);
        end
        
        % Get ts and tn
        ts = 0 - real(T);
        tn = p(:,ii) + imag(T) + .5*(Om(:,ii)+conj(Om(:,ii)));

        beta(ii,:) = AE_crack_solver(ts,tn,term(:,ii),A);

    end
    % Solver the tunnel
    T_func = @(z,ang,m_not_f,m_not_t) T_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,ang,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not_f,m_not_t);
    
    % Solve for the a and b coefficients (tunnel elements)
    for jj=1:length(zt)
        a(jj,:) = conj(1i*( Cauchy_integral_le( 1,N,m,max(m),...
            @(Z)z_from_bigz(Z,rt(jj),zt(jj)),...
            g,rho,rho_w,k,Omega_w_func,T_func,jj ) ));
        b(jj,:) = -conj( 1*Cauchy_integral_le( 0,N,m,max(m),...
            @(Z)z_from_bigz(Z,rt(jj),zt(jj)),...
            g,rho,rho_w,k,Omega_w_func,T_func,jj ) );
    end
    
    
    errorb = max(max(abs(b_old-b)./abs(b)));
    errorbeta = max(max(abs(beta_old-beta)./abs(beta)));
    errora = max(max(abs(a_old-a)./abs(a)));
    error = max([errora,errorb,errorbeta]);
    
    beta_old = beta;
    a_old = a;
    b_old = b;
    NIT = NIT + 1;
    
    % Display the error and iteration number
    disp_str = [' ',num2str(error,'%e'),'   ',num2str(NIT,'%i')];
    disp(disp_str)
end
end

