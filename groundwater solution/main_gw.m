function [Phi0,alpha,b,a,beta,c,t,Q] = main_gw (kb_in,z1b,z2b,ka_in,z1a,z2a,W_uni,rt,zt,z1gw,z2gw,z_ref,k)
%% Defintion of geometries

m = 10;
mfar = m*10;
m2 = m;
N = m;
Ne = m2*2;


bb = ones(1,length(z1b))./1;
ba = ones(1,length(z1a))./1;
kb = ones(1,length(z1b)).*k.*kb_in;
ka = ones(1,length(z1a)).*k./ka_in;
k_crack = [kb,ka];
b_crack = [bb,ba];

rho = 1000;
g = 10;
h = 0;

% Reference point
fi0 = imag(z_ref);
Phi0 = Phi_from_fi(fi0,k);

%% pre-calcualtions
beta = get_beta(m+mfar);
alpha = zeros(m,mfar);
for nn = 0:m-1
    alpha(nn+1,:) = get_alpha(beta,nn,mfar);
end

%% Iterator
[b,a,c,t,Q,Phi0] = iterator(W_uni,z1b,z2b,z1a,z2a,beta,alpha,m,m2,N,k,k_crack,b_crack,z_ref,Phi0,zt,rt,z1gw,z2gw,Ne);


%% Plot
xfrom = -1;
xto = -xfrom;
yfrom = -1;
yto = 1;
Nx = 100;
Ny = Nx;
lvs = 20;

xfrom2 = -500 - 50;
xto2 = 500 + 50;
yfrom2 = -530 - 50;
yto2 = 470 + 50;

load SF_COMSOL_chead_data
% load SF_COMSOL_chead2_data

for xx = 1:length(XC)
    z = complex(XC(xx),YC(xx));
    W = -W_total(z,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,rt,zt,t,Q,z1gw,z2gw,c,0,0,0,0);
    p(xx) = pressure( z,h,rho,g, @(z) Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,0,0,0,0), @(Phi) fi_from_Phi( Phi,k ) );
    SF(xx) = rho*g*W;
end
figure
hold on
plot(p)
th = linspace(0,2*pi,100);
create_figure_axis
plot(th,real(SF(1:100)),'-','Color',[.2 .5 0])
plot(th,SFXC(1:100),'--','Color',[.2 .5 0 ])
plot(th,-imag(SF(1:100)),'-','Color',[.5 0 .5])
plot(th,SFYC(1:100),'--','Color',[.5 0 .5])
legend('$s_{x,AEM}$','$s_{x,COMSOL}$','$s_{y,AEM}$','$s_{y,COMSOL}$','Location','southwest')
xticks([0, pi/2, pi, 3*pi/2, 2*pi])
xticklabels({'$0$','$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$' })
xlim([0 2*pi])
grid minor
ylabel('$s_j$')
create_figure_axis
plot(th,abs(SFXC(1:100)' - real(SF(1:100)))./(abs(SFXC(1:100)')),'Color',[.2 .5 0])
plot(th,abs(SFYC(1:100)' + imag(SF(1:100)))./(abs(SFYC(1:100)')),'Color',[.5 0 .5])
legend('$\Delta{}s_x$','$\Delta{}s_y$','Location','northwest')
xticks([0, pi/2, pi, 3*pi/2, 2*pi])
xticklabels({'$0$','$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$' })
xlim([0 2*pi])
grid minor
ylabel('$\Delta{}s_j$')

[Grid1,X1,Y1] = Contour_flownet(xfrom, xto, Nx, yfrom, yto, Ny, @(z)Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,0,0,0,0));
create_figure(600,600)
Contour_flow_net(X1,Y1,Grid1,lvs,1.5);
Plot_line(z1b,z2b,'black',1.5)
Plot_line(z1gw,z2gw,'black',1.5)
Plot_circle(zt,rt,'black',1.5)
axis([xfrom xto yfrom yto])
text(real(zt),imag(zt)+rt-3,'$\frac{\pi}{2}$','FontSize',16,'HorizontalAlignment','center' )
text(real(zt),imag(zt)-rt+3,'$\frac{3\pi}{2}$','FontSize',16,'HorizontalAlignment','center' )
text(real(zt)-rt+2,imag(zt),'$\pi$','FontSize',16,'HorizontalAlignment','center' )
text(real(zt)+rt-2,imag(zt),'$0$','FontSize',16,'HorizontalAlignment','center' )
[Grid12,X12,Y12] = Contour_flownet(xfrom2, xto2, Nx, yfrom2, yto2, Ny, @(z)Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,0,0,0,0));
create_figure(600,600)
Contour_flow_net(X12,Y12,Grid12,lvs,1.5);
Plot_line(z1b,z2b,'black',1.5)
Plot_line(z1gw,z2gw,'black',1.5)
Plot_circle(zt,rt,'black',1.5)
axis([xfrom2 xto2 yfrom2 yto2])
create_figure(600,600)
contourf(X12,Y12,fi_from_Phi(real(Grid12),k),lvs);
Plot_line(z1b,z2b,'black',1.5)
Plot_line(z1gw,z2gw,'black',1.5)
Plot_circle(zt,rt,'black',1.5)
axis([xfrom2 xto2 yfrom2 yto2])
create_figure(600,600)
contourf(X1,Y1,fi_from_Phi(real(Grid1),k),lvs);
Plot_line(z1b,z2b,'black',1.5)
Plot_line(z1gw,z2gw,'black',1.5)
Plot_circle(zt,rt,'black',1.5)
axis([xfrom xto yfrom yto])
