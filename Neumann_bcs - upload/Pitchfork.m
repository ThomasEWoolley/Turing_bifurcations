clear all
close all
clc
figure('Position',[0 0.05 1 2/3]);
subplot(2,3,[1:3])
hold on

load(['./hom1D/pt271.mat'])
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


for i=1:22
    load(['./1D',num2str(i),'/pt500.mat'])
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
    
end

xlabel('$L$')
ylabel('$\max(u)\, \& \,\min(u)$')
set(gca,'fontsize',15)
axis([0 1 0 5])

%%
for j=1:22
    Names=dir(['1D',num2str(j),'/pt*']);
    Indices=[];
    
    for k=1:length(Names)
    Indices(k)=str2num(Names(k).name(3:end-4));
    end
    Indices=sort(Indices);
    AllIndices(j)={Indices};
    iter=0;
    for i=Indices(2:end)
        iter=iter+1;
        load(['./1D',num2str(j),'/pt',num2str(i),'.mat'])
        Liter(iter,j)=p.u(end);
        
        StabValue(iter,j)=0;
        if(~ismember(p.branch(2,end-1),[-2 -1 1 2]) && p.branch(3,end-1)  ==0); StabValue(iter,j)=1;  end
        if( ismember(p.branch(2,end-1),[-2 -1 1 2]) && p.branch(3,end)==0); StabValue(iter,j)=1;  end
        
    end
end
xx=linspace(0,1,p.np);
%%


load(['./1D1/pt',num2str(AllIndices{1}(26+1)),'.mat']);
L=0.25;
subplot(2,3,[1:3])
xl=xline(L,'k','linewidth',3);
xl.Alpha=1;
subplot(2,3,4)
hold on
title(['$L=$ ',num2str(round(L,2))])
plot(xx,p.u(1:p.np),'linewidth',2,'color','k')
load(['./1D1/pt',num2str(AllIndices{1}(9+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',1,'color','k')
load(['./1D2/pt',num2str(AllIndices{2}(45+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',2,'color','k')
formatter

%%

load(['./1D2/pt',num2str(AllIndices{2}(43+1)),'.mat']);
L=p.u(end);
subplot(2,3,[1:3])
xl=xline(L,'r','linewidth',3);
xl.Alpha=1;
subplot(2,3,5)
hold on
title(['$L=$ ',num2str(round(L,2))])
plot(xx,p.u(1:p.np),'linewidth',2,'color','r')
load(['./1D2/pt',num2str(AllIndices{2}(25+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',1,'color','r')
load(['./1D3/pt',num2str(AllIndices{3}(2+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',1,'color','r')
formatter

%%
load(['./1D3/pt',num2str(AllIndices{3}(20+1)),'.mat']);
L=0.75;
subplot(2,3,[1:3])
xl=xline(L,'linewidth',3,'color',[0 0.7 0]);
xl.Alpha=1;
subplot(2,3,6)
hold on
title(['$L=$ ',num2str(round(L,2))])
plot(xx,p.u(1:p.np),'linewidth',2,'color',[0 0.7 0])
load(['./1D4/pt',num2str(AllIndices{4}(39+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',2,'color',[0 0.7 0])
load(['./1D5/pt',num2str(AllIndices{5}(2+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',1,'color',[0 0.7 0])
load(['./1D7/pt',num2str(AllIndices{7}(21+1)),'.mat']);
L=p.u(end);
plot(xx,p.u(1:p.np),'linewidth',1,'color',[0 0.7 0])
formatter
export_fig('C:\Users\smatw1\Dropbox\Turing_dirichlet\Pictures\Pitchfork_branch.png','-r300')
function formatter
axis([0 1 0 5])
set(gca,'fontsize',15)
ylabel('$u$')
xlabel('$x$')
end