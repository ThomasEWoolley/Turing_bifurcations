clear all
close all
clc
fs=15;
x=linspace(-.5,.5);

for n=0:1
a(n+1)=8*(-1)^(n+1)/(pi*(2*n+3)*(2*n+1)*(2*n-1));
c(:,n+1)=cos((2*n+1)*pi*x);
end
FT=c*a';

for n=0:10
a(n+1)=8*(-1)^(n+1)/(pi*(2*n+3)*(2*n+1)*(2*n-1));
c(:,n+1)=cos((2*n+1)*pi*x);
end
figure('position',[0 0 .5 1/3])
subplot(1,2,1)
plot(0:10,a,'d')
xlabel('$n$')
ylabel('Fourier coeficient')
set(gca,'fontsize',fs)
subplot(1,2,2)
plot(x,cos(pi*x).^2,'r')
hold on
plot(x,FT,'b--')
xlabel('$x$')
legend('$\cos(\pi x)^2$','Two term Fourier approximation','location','no')
set(gca,'fontsize',fs)

export_fig('../../Pictures/FT_approx.png','-r300')