function [] = PlotMain(x,xtu,xtd,phu,phd,d, fdom, Vnd, Vnb, f, f2)

N=length(Vnd(:,1));
Nb=length(Vnb(:,1));

figure('units','normalized','outerposition',[0 0 1 1])

subplot(4,3,[1,2 4,5])

%yyaxis left
ylu=abs(xtu);
yld=abs(xtd);
semilogy(fdom,ylu,'b')
hold on
semilogy(fdom,yld, 'r')
hold off
xlabel('Frequency domain [Hz]')
ylabel('Receptance [m/N]')
grid on
legend('Undamped', 'Damped')
%Modes excited

%find the 2 modes that colaborate most: the ones excited.
for i=1:N+Nb
[maxv(i), maxp(i)]=max(x(:,i));%Value and position of the maximium of every mode contribution
end
[temp,originalpos] = sort( maxv, 'descend' );
n = temp(1:2);
p =originalpos(1:2);%FRF

subplot(4,3,3)
  
    stem(Vnd(:,p(1)), 'r-o')
    hold on
    stem(Vnb(1:N,p(1)),'r-x')
    stem(Vnb(N+1:Nb,p(1)),'b-x')
    xlim([1 N+1])
    ylim([-0.55 0.55])
    label=num2str(d(p(1))^0.5);
    title(['Frequency: ',label,' [Hz]', '; Disks: -o ', '; Blades: -x ' ])
    hold off
 
 
subplot(4,3,6)
  
    stem(Vnd(:,p(2)), 'r-o')
    hold on
    stem(Vnb(1:N,p(2)),'r-x')
    stem(Vnb(N+1:Nb,p(2)),'b-x')
    xlim([1 N+1])
    ylim([-0.55 0.55])
    hold off 
    label=num2str(d(p(2))^0.5);
    title(['Frequency: ',label,' [Hz]', '; Disks: -o ', '; Blades: -x ' ])

% Phase only
subplot(4,3,[7,8])
plot(fdom,phu)
hold on
plot(fdom,phd)
hold off
legend('Phase')
xlabel('Frequency domain [Hz]')
ylabel('Phase [°]')
legend('Undamped', 'Damped')
grid on

%Force
subplot(4,3,12)
stem(f)
xlim([1 2*N])
xlabel('DOF')
title('Excitation Force: o ; Piezo Force: x')
hold on
stem(f2, '-x')

% Modes Contribution
subplot(4,3,9)
for i=1:N+Nb
semilogy(fdom,abs(x(:,i)))
hold on
end
title([num2str(2*N), ' modes contribution'])
grid on
hold off

% Transfer Function only
subplot(4,3,[10,11])
semilogy(fdom,ylu,'b')
hold on
semilogy(fdom, yld, 'r')
hold off
legend('Displacement')
xlabel('Frequency domain[Hz]')
ylabel('Displacement [m]')
%ylim([ (min(yl)) (max(yl))])
legend('Undamped', 'Damped')
grid on
end

