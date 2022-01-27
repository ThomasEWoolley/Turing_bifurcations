clear all
close all
clc

Names={'004','005','006','007','008','009','01','02'};
Cols=[0 0 0;1 0 0;0 1 0;0 0 1; 0.4660 0.6740 0.1880 ;0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560;0.8500 0.3250 0.0980];

figure('position',[0 0.1 .5 1/3])
hold on

for j=1:length(Names)
    n=getlatestfile(['./1D1',Names{j}]);
    try
        load(['./1D1',Names{j},'/',n])
        for i=1:length(p.branch(2,:))-1 % loop from first-point to last point
            lw=1;
            if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
            if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end

            plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'color',Cols(j,:),'Linewidth',lw);
        end
    end
    n=getlatestfile(['./1D2',Names{j}]);
    try
        load(['./1D2',Names{j},'/',n])
        for i=1:length(p.branch(2,:))-1 % loop from first-point to last point
            lw=1;
            if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
            if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end

            plot(p.branch(4,i:i+1),p.branch(8,i:i+1),'color',Cols(j,:),'Linewidth',lw);
        end
    end
end
%%
for j=1:length(Names)
    plo(j)=plot([0 1], -[1 1],'color',Cols(j,:));
end
xlabel('$L$')
ylabel('$\max(u)\, \&\, \min(u)$')
legend(plo,'$D_v=0.04$','$D_v=0.05$','$D_v=0.06$','$D_v=0.07$','$D_v=0.08$','$D_v=0.09$','$D_v=0.10$','$D_v=0.20$','location','neo')
axis([0 0.15 0 120])
set(gca,'fontsize',15)
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Side_peak_dependence.png','-r300')

%%
clear all
close all
clc
u1=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\L015_Dv0.05.txt');
u2=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\L015_Dv0.1.txt');
u3=load('C:\Users\smatw1\Dropbox\Turing_dirichlet\Maths\Comsol\Data\L015_Dv0.2.txt');
x=u1(:,1);
u1(:,1)=[];
u2(:,1)=[];
u3(:,1)=[];
figure('position',[0 0.1 .5 1/3])
subplot(1,3,1)
plot(x,u1(:,end))
title('$D_v=0.05$')
formatter
subplot(1,3,2)
plot(x,u2(:,end))
formatter
title('$D_v=0.1$')
subplot(1,3,3)
plot(x,u3(:,end))
formatter
title('$D_v=0.2$')


export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Side_peak_solutions.png','-r300')

function formatter
xlabel('$x$')
ylabel('$u$')
% legend('$D_v=0.05$','$D_v=0.01$','$D_v=0.02$','location','n')
axis([-0.5 0.5 0 80])
set(gca,'fontsize',15)
end