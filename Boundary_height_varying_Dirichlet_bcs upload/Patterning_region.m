clear all
close all
clc

[L,dv]=meshgrid(linspace(0,0.5,99),linspace(0,0.5));

lambda=(-(1./2).*pi.^2.*dv-3./2.*(L.^2)-(1./2000).*pi.^2+...
    (1./2000).*sqrt(1000000.*pi.^4.*dv.^2+...
    10000000.*L.^2.*pi.^2.*dv-2000.*pi.^4.*dv-23000000.*L.^4-10000.*L.^2.*pi.^2+pi.^4))./L.^2;
pcolor(L,dv,real(lambda))

shading interp
cbh=colorbar;
caxis([0 0.8]);
cbh.Ticks = 0:.2:.8; %Create 8 ticks from zero to 1
cbh.TickLabels = {'$\leq 0$','0.2','0.4','0.6','0.8'} ;
cbh.TickLabelInterpreter = 'latex';
xlabel('$L$')
ylabel('$D_v$')
set(gca,'fontsize',15)
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Patterning_region.png','-r300')

%%
clear all
close all
clc
dv=linspace(0,0.4,1e6);
m=3/10.*((sqrt(62500.*dv.^2-2500.*dv+1)-...
    250.*dv+5).*sqrt(2).*(62500.*sqrt(62500.*dv.^2-...
    2500.*dv+1).*dv.^2-15625000.*dv.^3-...
    2000.*dv.*sqrt(62500.*dv.^2-2500.*dv+1)+...
    812500.*dv.^2+3.*sqrt(62500.*dv.^2-2500.*dv+1)-7750.*dv+3)./...
    (sqrt(-5+1250.*dv-5.*sqrt(62500.*dv.^2-...
    2500.*dv+1)).*dv.*(250.*dv.*sqrt(62500.*dv.^2-...
    2500.*dv+1)-62500.*dv.^2-11.*sqrt(62500.*dv.^2-...
    2500.*dv+1)+4000.*dv-43).*(1-250.*dv+sqrt(62500.*dv.^2-2500.*dv+1))));
ind = find(abs(imag(m))<eps);
plot(dv(ind),m(ind))

xticks([0 0.05 0.1:0.1:0.5])
xlabel('$D_v$')
ylabel('$m$')
set(gca,'fontsize',15)
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Amplitude_steady_state_dependence_on_Dv.png','-r300')
