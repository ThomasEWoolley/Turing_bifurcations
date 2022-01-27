clear all
close all
clc
figure('Position',[0 0.05 1 0.8]);
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


load('./1D2/pt1845.mat')
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
axis([0 2 0 25])



u1=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\Subcritical_trajectory.txt');
x=u1(:,1);
u1(:,1)=[];

u2=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\Supercritical_trajectory.txt');
u2(:,1)=[];
u3=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\Supercritical_later_trajectory.txt');
u3(:,1)=[];

[m,n]=size(u1);
L1=0.1*(1+linspace(0,3000,n)/100);
L2=0.4*(1+linspace(0,3000,n)/100);
%%

p1=plot(L1,max(u1),'k--')
hold on
plot(L1,min(u1),'k--')
p2=plot(L1,max(u2),'r--')
plot(L1,min(u2),'r--')
p3=plot(L2,max(u3),'--','color',[0 0.7 0])
plot(L2,min(u3),'--','color',[0 0.7 0])
xlim([min(L1),2])
p4=plot([0 -1],[0 -1],'b','linewidth',1);
p5=plot([0 -1],[0 -1],'b','linewidth',3);
legend([p5,p4,p1,p2,p3],'Stable solutions','Unstable solutions','Subcritical branch','Supercritical branch', 'Second supercritical branch','location','e')


xlabel('$L$')
ylabel('$\max(u)\, \& \,\min(u)$')
set(gca,'fontsize',15)
xticks([ 0.1 0.2:0.2:2])
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Trajectory_bifurcations_Dirichlet.png','-r300')


%%

figure('position',[0 0 1 .3])
subplot(1,3,1)
plotfig(u1,L1,x)
title('Subcritical branch')
xticks([ 0.1 0.5 1 1.5 2])

subplot(1,3,2)
plotfig(u2,L1,x)
title('Supercritical branch')
xticks([ 0.1 0.5 1 1.5 2])

subplot(1,3,3)
plotfig(u3,L2,x)
title('Second supercritical branch')
xticks([0.4 1 1.5 2])

cbh=colorbar;
caxis([0 5]);
cbh.Ticks = 0:5; %Create 8 ticks from zero to 1
cbh.TickLabels = {'0','1','2','3','4','$>5$'} ;
cbh.TickLabelInterpreter = 'latex';
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Simulation_trajectories_Dirichlet.png','-r300')



function plotfig(u,L,x)
cmapvec=[0 5];
pcolor(L,x,abs(u))
shading interp
caxis(cmapvec)
xlim([min(L),2])
xlabel('$L$')
ylabel('$x$')
set(gca,'fontsize',15)
end