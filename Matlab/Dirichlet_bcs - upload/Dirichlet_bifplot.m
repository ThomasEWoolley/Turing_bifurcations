clear all
close all
clc
figure('Position',[0 0.05 1/2 1/3]);
hold on
pb=plot([0 1],[-1 -1],'md','linewidth',3);
p0=plot([0 1],[-1 -1],'g','linewidth',3);
p1=plot([0 1],[-1 -1],'color',[0 0.7 0],'linewidth',3);
p2=plot([0 1],[-1 -1],'k','linewidth',3);
p3=plot([0 1],[-1 -1],'b','linewidth',3);
p4=plot([0 1],[-1 -1],'r','linewidth',3);
%%
ln=linspace(-0.1,0,1e3);
lp=linspace(0,0.1,1e3);
ApproxBif1p=(13.92564203*lp/15.19200896);
ApproxBif1n=(13.92564203*ln/15.19200896);
plot(.1066440058+lp,2-30.26047934*ApproxBif1p,'g','linewidth',3)
plot(.1066440058+ln,2-30.26047934*ApproxBif1n,'g','linewidth',1);

epsilon=1;
l=linspace(-0.1,0.1,1e4);
ApproxBif21=-(0.1423206956e-1*(-268.0242070*epsilon*l+...
    15.19200896-sqrt((268.0242070*epsilon*l-15.19200896).^2+...
    140.5276998*epsilon*(-237.6901417*epsilon*l.^2+13.92564203*l))))/epsilon;
ApproxBif22=-(0.1423206956e-1*(-268.0242070*epsilon*l+...
    15.19200896+sqrt((268.0242070*epsilon*l-15.19200896).^2+...
    140.5276998*epsilon*(-237.6901417*epsilon*l.^2+13.92564203*l))))/epsilon;
ApproxBif2=[ApproxBif21,fliplr(ApproxBif22)];

StabA21=((-3.814539154*epsilon*l+...
    .2162137283).*sqrt(38434.92666*epsilon^2*l.^2-...
    6186.713865*epsilon*l+230.7971362)-...
    547.0085496*epsilon^2*l.^2+88.0497421*epsilon*l-3.284720897)/epsilon;
StabA22=((3.814539154*epsilon*l-...
    .2162137283).*sqrt(38434.92666*epsilon^2*l.^2-...
    6186.713865*epsilon*l+230.7971362)-...
    547.0085496*epsilon^2*l.^2+88.0497421*epsilon*l-3.284720897)/epsilon;
StabA2=[StabA21,fliplr(StabA22)];
ll=[l,fliplr(l)];
RealIndices=ApproxBif2==real(ApproxBif2);
ll(RealIndices==0)=nan;
ApproxBif2(isnan(ll))=[];
StabA2(isnan(ll))=[];
ll(isnan(ll))=[];



for i=1:length(ll)-1 % loop from first-point to last point
    lw=1;
    
    if(StabA2(i)<0); lw=3;  end
    plot(.1066440058+ll(i:i+1),2-30.26047934*ApproxBif2(i:i+1),'Linewidth',lw,'color',[0 0.7 0]);
    
end


%%
n1=1;
n2=50;
load('./hom1D/pt65.mat')
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'k','Linewidth',lw);
end




%%
n1=1;
n2=80;
load('./1D1/pt1048.mat')
% plot(p.branch(4,n1:n2),p.branch(8,n1:n2),'b','linewidth',3) %Plots branch max
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'b','Linewidth',lw);
end

%%
n1=1;
n2=80;
load('./1D2/pt1845.mat')
% plot(p.branch(4,n1:n2),p.branch(9,n1:n2),'r','linewidth',3) %Plots branch min
for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end
    
    plot(p.branch(4,i:i+1),p.branch(9,i:i+1),'r','Linewidth',lw);
end
%%
plot(.1066440058,2,'md','linewidth',3)
axis([0.04 0.2 0 8])
xlabel('$L$')
ylabel('$\max(u)\, \& \,\min(u)$')
set(gca,'fontsize',15)
h=legend([pb,p0,p1,p2,p3,p4],'Bifurcation point','Quadratic amplitude equation','Cubic amplitude equation','Zero steady state','Subcritical branch','Supercritical branch','location','eo');
set(h,'FontSize',15);
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Bifurcation_point_Dirichlet.png','-r300')
