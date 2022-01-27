clear all
close all
clc
figure('Position',[0 0.05 1 2/3]);
subplot(2,5,[1:5])
hold on



load('./1D1/pt1048.mat')
n1=1;
n2=length(p.branch(8,:));

% plot(p.branch(4,n1:n2),p.branch(8,n1:n2),'b','linewidth',3) %Plots branch max
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'b','Linewidth',lw);
    plot(p.branch(4,i:i+1),p.branch(9,i:i+1),'b','Linewidth',lw);
end
xlabel('$L$')
ylabel('$\max(u)\, \& \,\min(u)$')
set(gca,'fontsize',15)
axis([0 1 0 50])

%%
iter=0;
for i=0:10:1000
    iter=iter+1;
load(['./1D1/pt',num2str(i),'.mat'])
Liter(iter)=p.u(end);
end
%%



load('./1D1/pt460.mat')
L=p.u(end);
subplot(2,5,[1:5])
xl=xline(L,'k','linewidth',3);
xl.Alpha=1;
subplot(2,5,6)
title(['$L=$ ',num2str(round(L,2))])
p.pdeo.grid.plot(p.u(1:p.np),'color','k')
load('./1D1/pt110.mat')
L=p.u(end);
p.pdeo.grid.plot(p.u(1:p.np),'linewidth',1,'color','k')
formatter


load('./1D1/pt570.mat')
L=p.u(end);
subplot(2,5,[1:5])
xl=xline(L,'r','linewidth',3);
xl.Alpha=1;
subplot(2,5,7)
title(['$L=$ ',num2str(round(L,2))])
% plot(p.u(1:p.np))
p.pdeo.grid.plot(p.u(1:p.np),'color','r')
load('./1D1/pt20.mat')
L=p.u(end);
p.pdeo.grid.plot(p.u(1:p.np),'linewidth',1,'color','r')
formatter


load('./1D1/pt800.mat')
L=p.u(end);
subplot(2,5,[1:5])
xl=xline(L,'color',[0 0.7 0],'linewidth',3)
xl.Alpha=1;
subplot(2,5,8)
title(['$L=$ ',num2str(round(L,2))])
% plot(p.u(1:p.np))
p.pdeo.grid.plot(p.u(1:p.np),'color',[0 0.7 0])
formatter


load('./1D1/pt900.mat')
L=p.u(end);
subplot(2,5,[1:5])
xl=xline(L,'linewidth',3,'color',[0 0 0.4])
xl.Alpha=1;
subplot(2,5,9)
title(['$L=$ ',num2str(round(L,2))])
% plot(p.u(1:p.np))
p.pdeo.grid.plot(p.u(1:p.np),'color',[0 0 0.4])
formatter


load('./1D1/pt940.mat')
L=p.u(end);
subplot(2,5,[1:5])
xl=xline(L,'k','linewidth',3,'color',[1 0 1])
xl.Alpha=1;
subplot(2,5,10)
title(['$L=$ ',num2str(round(L,2))])
% plot(p.u(1:p.np))
p.pdeo.grid.plot(p.u(1:p.np),'color',[1 0 1])
formatter
shg
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Subcritical_branch.png','-r300')
function formatter
axis([-0.5 0.5 0 50])
set(gca,'fontsize',15)
ylabel('$u$')
xlabel('$x$')
end