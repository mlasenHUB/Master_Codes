figure
aux=MminB1M;
aux(aux==0)=nan;
for i=1:8
    plot(aux(:,i),'d', 'MarkerSize',10)
    hold on
    legendInfo{i} = ['Target mode: ', num2str(i)]; 
end
xlim([1 12])
ylim([-1.01 1.01])
xticks(1:1:12) 
ylabel('Force [N]')
xlabel('Piezo Position')
legend(legendInfo);
grid on
hold off
line([1,16],[1,1],'Color','red','LineStyle','--','LineWidth',2,'DisplayName','Force = 1')
line([1,16],[-1,-1],'Color','red','LineStyle','--','LineWidth',2,'DisplayName','Force = -1')
