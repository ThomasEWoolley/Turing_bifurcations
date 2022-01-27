ccc
m = 0;
x = linspace(0,1,200);
t = linspace(0,100);
L=0.07005;
% L=1;
sol = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,L),@pdex1ic,@pdex1bcd,x,t);
u = sol(:,:,1);
v = sol(:,:,2);
ux=linspace(0,max(u(:)));

% A surface plot is often a good way to study a solution.
subplot(1,3,1)
pcolor(x,t,u)
shading interp

xlabel('Distance x')
ylabel('Time t')

subplot(1,3,2)
plot(x,u(end,:))
hold on
plot(x,v(end,:))

subplot(1,3,3)
plot(sol(end,:,1),sol(end,:,2))
hold on
plot(ux,(2*ux-1)./ux.^2)
plot(ux,3./ux.^2)
axis([0 max(ux) 0 max(max(sol(:,:,2)))])
xlabel('$u$')
ylabel('$v$')


% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx,L)
c = [1;1];
f = [1/1000;1/10].*DuDx/L^2;
s = [1-2*u(1)+u(1)^2*u(2);3-u(1)^2*u(2)];
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = [2*(1+x*(1-x)/0.25);0.75*(1+x*(x-1)/0.25)];
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bcd(xl,ul,xr,ur,t)
pl = ul-[2;0.75];
ql = [0;0];
pr = ur-[2;0.75];
qr = [0;0];
end
function [pl,ql,pr,qr] = pdex1bcn(xl,ul,xr,ur,t)
pl = [0;0];
ql = [1;1];
pr = [0;0];
qr = [1;1];
end