function [ ] = PlotTF( x,xt, f2, pp, fdom, r )
%Plot transfer functions and its mode contribution
%   Detailed explanation goes here
N=length(f2);
l=length(fdom);
figure('units','normalized','outerposition',[0 0 1 1])

for i=1:N
semilogy(fdom,abs(x(:,i)))
[maxv(i), maxp(i)]=max(x(:,i));%Value and position of the maximium of every mode contribution
hold on
legendInfo{i} = ['Mode: ', num2str(i)];  

end
grid on

semilogy(fdom,abs(xt),'-x')
legendInfo{N+1} = ('Total FRF'); 
legend(legendInfo);
 
 hold off
 xlabel('Frequency Domain [ Hz]')
 ylabel('FRF [m/N]')
 title(['Piezo pos: ',num2str(pp),'; No Damping in the model, Response measured at point ', num2str(r)])


end

