function [] = PlotTFpoints( c, xr, xtr,phr, fdom, h)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
n=length(xtr(1,:));
text=h;
bb=n/c;
%{
figure('units','normalized','outerposition',[0 0 1 1])
suptitle(['All response functions' ])

for i=1:n

subplot(bb,c,i)       
semilogy(fdom,abs(xtr(:,i)), 'r-x')
hold on
    for j=1:2*N
        semilogy(fdom,abs(xr{i}))
    end
hold off
label=(['Response point: ', num2str(i)]);
title(label)
xlabel('Frequency Domain')
ylabel('Displacement')
end
%}
figure('units','normalized','outerposition',[0 0 1 1])
suptitle(['FRF and Phases: ', text]);

subplot(2,1,1)
title('FRF for all points')

for i=1:n
   
semilogy(fdom,abs(xtr(:,i)))
hold on
legendInfo{i} = ['Response at point: ', num2str(i)];  
end
hold off
legend(legendInfo);
xlabel('Frequency Domain [Hz]')
ylabel('Amplitude')


subplot(2,1,2)
title('Phase for all points')
ylim([-200 200])

for i=1:n
   
plot(fdom,phr(:,i))
hold on
legendInfo{i} = ['Phase at point: ', num2str(i)];  
end
hold off
legend(legendInfo);
xlabel('Frequency Domain [Hz]')
ylabel('Phase')

end

