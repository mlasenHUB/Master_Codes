
%  3 M Machine System Dynamics   Tutorial 1 qu 5   DAR 1.10.00 
%frf regeneration   2  resonances uses complex variables
f=[1:500]; % freq range in Hz 1 Hz increments
m1=1; % mode 1  mass    in kg
f1= 100;  % mode 1 freq in hz
c1=0.1 ;  %  mode 1 viscous damping ratio 
m2= 1 ;  % mode 2 mass in kg
f2 = 200 ; %mode 2 freq in hz
c2=0.1 ;  %  mode 2 viscous damping ratio
w=2*pi*f; % rad /sec
w1=2*pi*f1; % res freq 1 , rad/sec
w2=2*pi*f2; % res freq 2 , rad/sec
k1=m1*w1*w1; % hence stiffness k1 in N/m

k2=m2*w2*w2; % hence stiffness k2 in N/m
x=w/w1;    %freq ratio     Mode 1 
a=(1-x.*x); % real denominator
b=c1.*x ;  % imag denominator
denom = complex(a,b);
h1=1/k1./denom;
rrecep=real(h1); 
irecep =imag(h1);
mag1=abs(h1);
phase =atan2(irecep,rrecep); % correct to fourth quadrant
phase=phase.*180/pi;  %convert to degrees
x2=w/w2;    %freq ratio     Mode 2 
a2=(1-x2.*x2); % real denominator
b2=c2.*x2 ;  % imag denominator
denom2=complex(a2,b2);
h2=1/k2./denom2;
rrecep2=real(h2); 
irecep2 =imag(h2);
mag2=abs(h2);
phase2 =atan2(irecep2,rrecep2); % correct to fourth quadrant
phase2=phase2.*180/pi;  %convert to degrees
% sum modes
rrecept=rrecep+rrecep2;
irecept=irecep+irecep2;
magt=sqrt(rrecept.*rrecept+irecept.*irecept);
phaset =atan2(irecept,rrecept); % correct to fourth quadrant
disp(' Two degrees of freedom model > Frequency Response Function' )
disp('  ')
disp('  Mode 1  mass 1 kg resonance  100 Hz,  0.1 viscous damping ratio. ')
disp('  Mode 2  mass 1 kg resonance  200 Hz,  0.1 viscous damping ratio. ')
disp('  ')
disp('  Frequency resolution:   1 Hz from 1 - 500 Hz.')
disp('   ')
disp(' Press return key to display FRF data in various forms')
pause
figure(1)
clf   % clear graphics
plot(f,mag1,'r',f,mag2,'g')
xlabel('frequency Hz.')
ylabel('Receptance Mag.')
title('Linear Receptance.Each mode ')
legend(' Mode 1','Mode 2') 
grid on
pause
figure(2)
clf
semilogy(f,mag1,'r',f,mag2,'g')
xlabel('frequency Hz.')
ylabel('Receptance Log Mag.')
title('Log Receptance. Each mode')
legend(' Mode 1','Mode 2') 
grid on
pause
figure(3)
clf
semilogy(f,mag1,'r',f,mag2,'g',f,magt,'b')
xlabel('frequency Hz.')
ylabel('Receptance Log Mag.')
title('Log Receptance. Sum')
legend('mode 1','mode2','sum')
% grid on
pause
figure(4)
clf
subplot(2,1,1)
semilogy(f,magt,'r')
xlabel('frequency Hz.')
ylabel('Receptance Log Mag.')
title('Log Receptance.')
grid on
subplot(2,1,2)
plot(f,phaset,'g')
xlabel('frequency Hz.')
ylabel('Phase degrees')
title('Phase angle')
grid on
pause
figure(5)
clf
plot (rrecept,irecept,'r')
xlabel(' Real Part')
ylabel ('imag part')
title(' Nyquist Plot of SDOF')
grid on 
pause
figure(6)
clf
plot3(f,rrecept,irecept,'r')
grid on