% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   An anlytical Element Model for Seepage Forces
%
%   By: Erik Toller
%   2023-10-13
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
close all
clearvars
clc
figure_font

addpath('groundwater solution')

%% INPUT
gamma_w = 10;
nu = .3;
G = 2*10^4;
kappa = 3 - 4*nu;
sigma_11 = 0;
g = 10;
rho = 2500;
rho_w = 1000;
gamma_w = rho_w*g;
h = 0;

%% Groundwater Data
% Fractures
z1b = -.5 - .5i;
z2b = .5+.5i;
kb_in = 1000;
% z1b = [];
% z2b = [];

ka_in = 10^1000;
z1a = []
z2a = [];

% Tunnel
rt = 7;
zt = [];

% Groundwater surface
z1gw = [];
z2gw = [];
z_ref = -500-530i-470i;

% Global properties
k = 1e-7*10^(7*0);
W_uni  = -1i*1i*k*1e-3;


%% GET GROUNDWATER DATA
[Phi0,alpha_gw,b_gw,a_gw,beta_gw,c_gw,t,Q] = main_gw(kb_in,z1b,z2b,ka_in,z1a,z2a,W_uni,rt,zt,z1gw,z2gw,z_ref,k);

%% LINEAR ELASTICITY
% FRACTURES
if length(z1b) == 4
    z1 = [z2b(1),z1b(3)];
    z2 = [z2b(2),z1b(4)];
elseif length(z1b) == 1
    z1 = [z1b(1)];
    z2 = [z2b(1)];
else
    z1 = [];
    z2 = [];
end


%% Pre-calcullation
L = zeros(1,length(z1));
mu = zeros(1,length(z1));
for ii = 1:length(z1)
    L(ii) = sqrt((z2(ii)-z1(ii))*conj(z2(ii)-z1(ii)));
    mu(ii) = angle(z2(ii)-z1(ii));
end

Omega_w_func = @(z) Omega_total(z,Phi0,W_uni,z1b,z2b,b_gw,z1a,z2a,a_gw,beta_gw,alpha_gw,length(b_gw),z_ref,rt,zt,t,Q,z1gw,z2gw,c_gw,0,0,0,0);
fi_ref = real(Omega_w_func(z_ref))/k;

%% SOLVE THE CEOFFICENTS
m = 10;
N = m*2;
cond = 1e-3;
[beta,a,b,theta,p1,z] = solv_coef(cond,z1,z2,L,m,N,mu,rt,zt,sigma_11,g,rho,rho_w,k,h,Omega_w_func,0);

%% CHECK THE SOLUTION
N_check = N;
theta_check = linspace(1e-7,2*pi-1e-7,N_check);
for ii = 1:N_check
    for jj = 1:length(z1)
        ang = mu(jj);
        z = z_from_chi(exp(1i*theta_check(ii)),z1(jj),z2(jj));
        p(ii) = pressure( z,h,rho_w,g,Omega_w_func, @(Phi) fi_from_Phi( Phi,k ) );
        T_Frac(jj,ii) = 1i.*p(ii) -  T_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,ang,sigma_11,g,rho,rho_w,k,Omega_w_func,0,0);
    end
    for jj = 1:length(zt)
        ang = theta_check(ii) + pi/2;
        z = z_from_bigz(exp(1i*theta_check(ii)),rt(jj),zt(jj));
        T_Tunnel(jj,ii) = T_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,ang,sigma_11,g,rho,rho_w,k,Omega_w_func,0,0);
    end
end

if ~isempty(z1b)
create_figure_axis(900,600)
xticks([0 pi/4 pi/2 3*pi/4 pi 5*pi/4 6*pi/4 7*pi/4 2*pi])
xticklabels({'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
grid on
xlim([0 2*pi])
plot(theta_check,real(T_Frac))
plot(theta_check,-imag(T_Frac))
legend('$\Re T$','$\Im T$')
end
if ~isempty(zt)
create_figure_axis(900,600)
xticks([0 pi/4 pi/2 3*pi/4 pi 5*pi/4 6*pi/4 7*pi/4 2*pi])
xticklabels({'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
grid on
xlim([0 2*pi])
plot(theta_check,real(T_Tunnel))
plot(theta_check,-imag(T_Tunnel))
legend('$\Re T$','$\Im T$')
end
%% PLOT
Nx = 80;
Ny = Nx;
lvs = 10;
Nx_quiv = 20;
Ny_quiv = Nx_quiv;
Ntraj = 200;
lvs_traj = 20;
Ntrajw = 200;
lvs_trajw = Nx_quiv;


dim_plot = [-1,1,-1,1];
ds = 2/(Ntraj)*1.5;

% FRACTURE AND TUNNEL CLOSE
xfrom2 = -.5-.03;
xto2 = -.5+.03;
yfrom2 = -.5-.03;
yto2 = -.5+.03;
dim_plot_close = [xfrom2,xto2,yfrom2,yto2];

[Grid_1,Grid_2,Grid_t,X,Y] = Contour_major_minor(dim_plot(1), dim_plot(2), Nx, dim_plot(3), dim_plot(4), Ny, @(z) major_minor_stress(z,@(z) tau_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,sigma_11,g,rho,rho_w,k,Omega_w_func,0,0)));


create_figure()
contour(X, Y, Grid_1,lvs,'blue');
Plot_line(z1,z2,'black',1)
Plot_circle(zt,rt,'black',1)
% legend('$\sigma_{1}$')
axis(dim_plot)
%PLOT BOX
zz1 = [cv(xfrom2,yfrom2),cv(xfrom2,yto2),cv(xto2,yto2),cv(xto2,yfrom2)];
zz2 = [cv(xfrom2,yto2),cv(xto2,yto2),cv(xto2,yfrom2),cv(xfrom2,yfrom2)];
Plot_line(zz1,zz2,'black :',1.5)
% 
% 
create_figure()
contour(X, Y, Grid_2,lvs,'red');
Plot_line(z1,z2,'black',1)
Plot_circle(zt,rt,'black',1)
% legend('$\sigma_{2}$')
axis(dim_plot)
%PLOT BOX
zz1 = [cv(xfrom2,yfrom2),cv(xfrom2,yto2),cv(xto2,yto2),cv(xto2,yfrom2)];
zz2 = [cv(xfrom2,yto2),cv(xto2,yto2),cv(xto2,yfrom2),cv(xfrom2,yfrom2)];
Plot_line(zz1,zz2,'black',1)

[Grid_1_close,Grid_2_close,Grid_t_close,XX,YY] = Contour_major_minor(dim_plot_close(1), dim_plot_close(2), Nx, dim_plot_close(3), dim_plot_close(4), Ny, @(z) major_minor_stress(z,@(z) tau_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,sigma_11,g,rho,rho_w,k,Omega_w_func,0,0)));


create_figure()
contour(XX, YY, Grid_1_close,lvs,'blue');
Plot_line(z1,z2,'black',1)
Plot_circle(zt,rt,'black',1)
% legend('$\sigma_{1}$')
axis(dim_plot_close)

create_figure()
contour(XX, YY, Grid_2_close,lvs,'red');
Plot_line(z1,z2,'black',1)
Plot_circle(zt,rt,'black',1)
% legend('$\sigma_{2}$')
axis(dim_plot_close)

