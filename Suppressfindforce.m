function [Mmin2, Mnom2] = Suppressfindforce( f, V )
%Mmin2: What is the best force magnitude, to place at any position, to suppress a certain mode
a=f*V;
n=length(a);
Mmin2=zeros(n,n);
for k = 1:n%Mode Shape
    for ii =1:n%Forcing pos
        Mmin2(ii, k)= -a(k)/V(ii, k);
    end    
end

Mmin2(~isfinite(Mmin2))=NaN;%Neglect infinite values of forcing
Mmin2(isnan(Mmin2))=0;
Mnom2=Mmin2./max(abs(Mmin2(:)));
%{
figure('units','normalized','outerposition',[0 0 1 1])

for i=1:2*N
    plot(Mnom2(:,i), '--o')
    hold on
    legendInfo{i} = ['mode shape: ' num2str(i)];    
end

xlabel('Forcing position')
ylabel('Normalised Force needed to suppress')
title(['Force and position needeed to suppress each of the ', num2str(2*N), ' mode shpaes'])

%Standard deviation of the needeed forces per position, we want to minimize
%the necessary force to suppress all the modes
hold on
for i=1:2*N
stdv(i)=std(Mnom2(i,:));
end
plot(stdv, '--r*')
hold off
legendInfo{2*N+1} = ['Standard deviation of the forces'];
legend(legendInfo)
%}

end

