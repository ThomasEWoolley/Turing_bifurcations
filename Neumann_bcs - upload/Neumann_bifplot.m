clear all
close all
clc
figure('Position',[0 0.05 1/2 1/3]);
hold on
pb=plot([0 1],[-1 -1],'md','linewidth',3);
p0=plot([0 1],[-1 -1],'k','linewidth',3);
p1=plot([0 1],[-1 -1],'color',[0 0.7 0],'linewidth',3);

p2=plot([0 1],[-1 -1],'b','linewidth',3);

%%

lp=linspace(0,0.1,1e3);
ApproxBif1p=sqrt(13.92547685*lp/292.523001);
ApproxBif1n=-ApproxBif1p;
plot(.1066440058+lp,2-30.26047934*ApproxBif1p,'color',[0 0.7 0]','linewidth',3)
plot(.1066440058+lp,2+30.26047934*ApproxBif1p,'color',[0 0.7 0],'linewidth',3)



%%
n1=1;
n2=50;
load('./hom1D/pt271.mat')
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'k','Linewidth',lw);
    
end




%%
n1=1;
n2=80;
load('./1D1/pt500.mat')
% plot(p.branch(4,n1:n2),p.branch(8,n1:n2),'b','linewidth',3) %Plots branch max
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'b','Linewidth',lw);
    plot(p.branch(4,i:i+1),p.branch(9,i:i+1),'b','Linewidth',lw);
end


%%
plot(.1066440058,2,'md','linewidth',3)
axis([0.04 0.2 0 5])
xlabel('$L$')
ylabel('$\max(u)\, \& \,\min(u)$')
set(gca,'fontsize',15)
h=legend([pb,p1,p0,p2],'Bifurcation point','Amplitude equation','Zero steady state','Supercritical pitchfork','location','eo');
set(h,'FontSize',15);
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Bifurcation_point_Neumann.png','-r300')
% 