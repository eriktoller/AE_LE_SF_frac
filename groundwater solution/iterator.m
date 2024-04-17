function [b,a,c,t,Q,Phi0] = iterator(W_uni,z1b,z2b,z1a,z2a,beta,alpha,m,m2,N,k,k_crack,b_crack,z_ref,Phi0,zt,rt,z1gw,z2gw,Ne)
% Set the initial values
error = 1;
b = zeros(length(z1b),m);
b_old = b;
a = zeros(length(z1a),m);
a_old = a;
Tb = zeros(length(z1b),1);
t = zeros(length(zt),m2+1);
t_old = ones(length(zt),m2+1);
c = zeros(length(z1gw),m2+1);
c_old = ones(length(z1gw),m2+1);
Q_old=zeros(1,length(zt)+length(z1gw));
if isempty(z1b)
    b = [];
    b_old = [];
end
if isempty(zt)
    t = [];
    t_old = [];
end
if isempty(z1gw)
    c = [];
    c_old = [];
end
NIT = 0;

% Calcualte the Am martix
[AM] = get_AMQ(N,rt,zt,z1gw,z2gw,z_ref);
Phi_tunnel = Phi_from_fi(imag(zt),k);
Phi_gw_level = Phi_from_fi(imag(imag((z1gw+z2gw)./2)),k);
Phi_gw_level = [Phi_from_fi(1000,k),Phi_from_fi(0,k)];

% Calcualte the evaluation points
deltheta = 2*pi/N;
theta = [linspace(deltheta/2,pi-deltheta,N/2),linspace(pi+deltheta/2,2*pi-deltheta,N/2)];
theta = linspace(0+deltheta,pi-deltheta,N);
chi = exp(1i*theta);

% Get teh fraction for jump in hydraulic conductivity
for bb = 1:length(z1b)
    Tb(bb) = k/(k_crack(bb)*b_crack(bb));
end
for aa = 1:length(z1a)
    Ta(aa) = k_crack(length(z1b)+aa)/(k*b_crack(length(z1b)+aa));
end
% Calcualte the A martices
for kk = 1:length(z1b)
    % Get the length and the Theta constant
    L = sqrt((z2b(kk)-z1b(kk))*conj(z2b(kk)-z1b(kk)));
    Thetaf = 4.*chi.^2./( L.*(chi.^2-1) )./(2*pi);
    % Get the values for A
    AA = zeros(N+2,m);
    for ii = 1:N
        for jj = 1:m
            nn = jj-1;
            a1 = F_prim_func(chi(ii),beta,jj-1)*Thetaf(ii);
            AA(ii,jj) = real(a1) - Tb(kk)*cos(nn*theta(ii));
        end
    end
    AA(N+1,:) = cos((0:m-1)*pi);
    AA(N+2,:) = -1;
    % Save the values of A to a cell
    Ab{kk} = AA;
end
for kk = 1:length(z1a)
    % Get the length and the Theta constant
    L = sqrt((z2a(kk)-z1a(kk))*conj(z2a(kk)-z1a(kk)));
    Thetaf = 4.*chi.^2./( L.*(chi.^2-1) )./(2i*pi);
    % Get the values for A
    AA = zeros(N+2-2,m);
    for ii = 1:N
        for jj = 1:m
            nn = jj-1;
            a1 = F_prim_func(chi(ii),beta,jj-1)*Thetaf(ii);
            AA(ii,jj) = imag(a1) + Ta(kk)*cos(nn*theta(ii));
        end
    end
    % Save the values of A to a cell
    Aa{kk} = AA;
end

if ~isempty(z1b)
    posb  = pos_bulider( z1b,z2b );
end
if ~isempty(z1a)
    posa  = pos_bulider( z1a,z2a );
end

str = '     Error    Iteration';
disp(str)

while error > 1e-03 && NIT < 200
    % Get the discharges
    [Q] = solve_Q( z_ref,Phi0,Phi_tunnel,Phi_gw_level,AM,Ne,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,rt,zt,t,z1gw,z2gw,c );
    C = Q(end);
    
    % Solve for the b coefficients (draining fractures)
    for ii = 1:length(z1b)
        b(ii,:) = solve_bcoef(theta,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,N,Ab{ii},posb{ii},rt,zt,t,Q,z1gw,z2gw,c,ii);
    end
    
    % Solve for the a coefficients (blocking fractures)
    for ii = 1:length(z1a)
        a(ii,:) = solve_acoef(theta,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,N,Aa{ii},posa{ii},rt,zt,t,Q,z1gw,z2gw,c,ii);
    end
    
    % Solve for the d coefficients (tunnel elements)
    for ii=1:length(zt)
        t(ii,:) = -conj(Cauchy_integral_pressure( N,m2,max(m2),rt(ii),zt(ii),k,...
            @(Z)z_from_bigz(Z,rt(ii),zt(ii)),...
            @(z)Omega_total(z,C,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,0,0,ii,0) ));
    end
    
    % Solve for constant gw level
    for ii=1:length(z1gw)
        c(ii,:) = -conj(Cauchy_integral( N,m2,max(m2),...
            @(chi)z_from_chi(chi,z1gw(ii),z2gw(ii)),...
            @(z)Omega_total(z,C,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,0,0,0,ii) ));
    end
    
    % Get the maximum error
    errorb = max(max(abs(b_old-b)./abs(b)));
    errora = max(max(abs(a_old-a)./abs(a)));
    errort = max(max(abs(t_old-t)./abs(t)));
    errorc = max(max(abs(c_old-c)./abs(c)));
    errorQ = max(max(abs(Q_old-Q)./abs(Q)));
    if max(max(abs(a_old-a)))<1e-10
        errora = 0;
    end
    if max(max(abs(c_old-c)))<1e-10
        errorc = 0;
    end
    error = max([errorb,errora,errort,errorc,errorQ]);
    
    % Add to the iteration counter
    NIT = NIT + 1;
    
    % Reassing solution for a and b to storage
    b_old = b;
    t_old = t;
    c_old = c;
    a_old = a;
    Q_old = Q;
    
    % Display the error and iteration number
    disp_str = [' ',num2str(error,'%e'),'   ',num2str(NIT,'%i')];
    disp(disp_str)
    
end
Phi0 = C;

end

