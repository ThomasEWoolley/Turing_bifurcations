clear all
close all
clc

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

%%
close all
hold on
which=7;
plot(Liter(:,which).*StabValue(:,which),'bd')
plot(Liter(:,which).*~StabValue(:,which),'ks')