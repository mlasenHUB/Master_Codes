function [  ] = PlotOpt(c,Ft, Pt , Mf,h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
text=h;
[r,col,deep]=size(Ft);
mf=zeros(1, r);%vector of maximum force among each amount of number of piezos 1...N
minf=zeros(1,r);
Ptf=Pt(r,:);
Ptfv=zeros(r, r); %for visualisation
for i =1 :r
    if isnan(Ptf(i))
        break
    end
    
    Ptfv(Ptf(i),i)=1;
    if i>1
        Ptfv(:,i)=Ptfv(:,i)+Ptfv(:,i-1);
    end
end

for i=1:r
    mf(i)=max(Mf(3*i-2,:));
    minf(i)=max(Mf(3*i-1,:));
end


figure('units','normalized','outerposition',[0 0 1 1])
suptitle(['Optimization for 1 piezo per mode suppress - ',text]);

colormap(flipud(gray))
imagesc(Ptfv)
colorbar
xlabel('Amount of piezos')
yyaxis left
ylabel('Piezo Position')


yyaxis right
plot(mf,'-o')
ylabel('Max and Min forces among the piezos')
set(gca,'xtick',linspace(1,r,r))
hold on
plot(minf,'-x')
grid on

figure('units','normalized','outerposition',[0 0 1 1])
suptitle(['Total force per mode for all possible piezos - 1 piezo per mode suppression - ',text]);

bb=r/c;
for i=1:r
   

subplot(bb,c,i)
    for j=1:i
        
        legendInfo{j} = ['Piezo ', num2str(j)];
        stem(Ft(j,:,i))
        legend(legendInfo);
        set(gca,'xtick',linspace(1,r,r))
        xlabel('Modes')
               
        hold on
    end
    
end
end

