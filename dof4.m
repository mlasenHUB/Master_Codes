clear all
close all

%%CONTENTS
%(1)Definition of K and M for 24 blades+disks (48 DOF)
%(2)Solution of eigen problem
%(3)Plotting of mode shapes nh

%%
%K-STIFFNESS MATRIX
N=4;
kdr=1000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kb=1000*ones(N,1);%stiffness of the disk-blade spring of disk-blade 'i'
kg=1000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'
z=zeros(N,1);

d= kdl+kdr+kb+kg;%elements of K from the equations of motion of the disks
d= [d;kb];%concatenate elements from equations of motion of the blades

KD= diag(d);%create basic Stiffness matrix diagonal, without coupling elements 
Kr=diag([(-kdr(1:N-1));z],1);%to the right of the diagonal
Krc=diag((-kb), N);%to the right corner of the diagonal
K=Kr+Krc;
K(1,N)=-kdr(N);%term due to the left of the first disk
K=K+K.';%Symmetric
K=K+KD;%STIFFNESS MATRIX 

%%
%M-Mass Matrix
md=ones(N,1);%masses of the disks
mb=ones(N,1);%masses of the blades
di=[md;mb];
M=diag(di);%MASS MATRIX


%%
%EIGEN PROBLEM
[V,D] = eig(K,M);%V:eigen vectors, D:eigen values.

%%
%MASS NORMALISATION
MR=V'*(M*V);%modal mass calculation
KR=V'*(K*V);%modal stiffness calculation
Vn=V*(MR^(0.5));%Normalised eigenvectors
d=diag(D);%just the (squared) frequencies

for i=1:2*N
    len(i)=norm(Vn(:,i));
end
%%
%MODE SHAPES
Vnd=Vn(1:N,:);%Disks Modal Shapes
Vnb=Vn(N+1:2*N,:);%Blades Modal Shapes

%%
%PLOTTING 48 Mode Shapes
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
subplot(4,4,i)       
stem(Vnd(:,i), 'r')
hold on
stem(Vnb(:,i),'x')
xlim([1 N])
%ylim([-0.3 0.3])
label=num2str(d(i)^0.5);
title(label)
%set(get(gca,'title'),'Position',[5.5 0.2 0.5])
hold off
end
legend('Disk', 'Blades')
suptitle('All Mode Shapes')

%%
%TRANSFER FUNCTIONS

e=1;% Excitation Point
r=2;% Response Point/Drive Point
fdom=[1:0.01:1.5*(max(D(:))^0.5)];%Frequency domain
k=length(d);

for i=1:k
    A(i)=Vn(e,i)*Vn(r,i);%Modal Constant
end
for j=1:k
    for i=1:length(fdom)
        H(i,j)=A(j)/((d(j)-fdom(i)^2));%Mode contribution
    end
end


%%
%TRANSFER FUNCTIONS FOR EVERY r AND  e IN THE SYSTEM
%{
for r=1:2*N%row: excitation points 
    for c=1:2*N%column: response points
        for i=1:k
            A(i)=Vn(r,i)*Vn(c,i);
        end
        for j=1:k
            for i=1:length(fdom)
                H(i,j)=A(j)/(d(j)-fdom(i)*fdom(i));
            end
        end
        H=abs(H);
        for i=1:length(H(:,1))
            Ht(i)=sum(H(i,:));%sum all the mode's contribution
        end
        Ha{r,c}=Ht;
    end
end
%}
%%
%PLOTTING TRANSFER FUNCTIONS

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(H(:,i)))
hold on
end
title('All 48 modes contribution')
grid on
hold off

for i=1:length(fdom)
    Ht(i)=sum(H(i,:));
end

figure('units','normalized','outerposition',[0 0 1 1])
%figure
semilogy(fdom,abs(Ht))
grid on
title('Total transference function, No Damping')
%{
descr = {'Excitation point' num2str(e);
    'Response point' num2str(r)};
ax2 = axes('Position',[.3 .1 .6 .8]);
axes(ax2) % sets ax1 to current axes
text(.025,0.6,descr)
%}
%%
%PRINTING ANY TRANSFER FUNCTION
%{
figure('units','normalized','outerposition',[0 0 1 1])
semilogy(fdom,H{e,r})
grid on
title('Total transference function')
%}
%%
%(1.a)Damping Models: SIMPLE PROPORTIONAL DAMPING C=bK.

b=0.0005;%proportional coefficient
kr=diag(KR);%Diagonal of modal stiffnes matrix.
mr=diag(MR);%Diagonal of modal mass matrix.
cr=b.*kr;%Diagonal of modal damping matrix.
for i=1:k
    Aspd(i)=V(e,i)*V(r,i);%Not normalized eigenvectors, SPD:Simple Prop. Damp.
end
for j=1:k
    for i=1:length(fdom)
        den=complex((kr(j)-mr(j)*fdom(i)^2), fdom(i)*cr(j));%denominator.
        Hspd(i,j)=Aspd(j)/den;%Mode contribution SPD:Simple Prop. Damp.
    end
end
Hspda=abs(Hspd);%Absolute value

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(Hspda(:,i)))
%plot(fdom,Hspdr(:,i))
hold on
end
title('All 48 modes contribution with Simple Proportional Damping')
grid on
hold off
for i=1:length(fdom)
    Htspdc(i)=sum((Hspd(i,:)));%Total Complex transfer function.
end
for i=1:length(fdom)
    Htspda(i)=abs(Htspdc(i));%Absolute part of the total transfer function.
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
semilogy(fdom,abs(Htspda))
grid on

hold off
subplot(1,2,2)
plot(real(Htspdc), imag(Htspdc))
grid on
suptitle('Total transference function Simple Proportional Damping')


%%
%(1.b)Damping Models: PROPORTIONAL DAMPING, HYSTERIC C=H=b2K+gM
b2=0.0005;%proportional coefficient- Stiffness
g=0.05;%proportional coefficient- Mass
nr=b2+g./d;%damping loss factor
cr=b*kr;%Diagonal of modal damping matrix.
for i=1:k
    Ah(i)=V(e,i)*V(r,i);%Not normalized eigenvectors, h:hysteric
end
for j=1:k
    for i=1:length(fdom)
        denhd=complex((kr(j)-mr(j)*fdom(i)^2), nr(j)*kr(j));%denominator.
        Hh(i,j)=Ah(j)/denhd;%Mode contribution SPD:Simple Prop. Damp.
    end
end
Hha=abs(Hh);%Absolute part

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(Hha(:,i)))
hold on
end
title('All 48 modes contribution with Hysteric Damping')
grid on
hold off
for i=1:length(fdom)
    Hhtc(i)=sum((Hh(i,:)));%Complex total transfer function.
end

for i=1:length(fdom)
    Hhta(i)=abs(Hhtc(i));%Absolute part of the total transfer function.
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
semilogy(fdom,abs(Hhta))
grid on

hold off
subplot(1,2,2)
plot(real(Hhtc), imag(Hhtc))
grid on
suptitle('Total transference function Hysteric Damping')
%%
%(1.c)Damping Models: NON- PROPORTIONAL DAMPING ??? PENDING

%%
%(2.a)Forcing Models: MULTIPLE FORCES, SAME MAGNITUDE. and PHASE.(Hysteric
%Damping)
%Creating an excitation in any/various DOF.
l=1;
for i=1:2*N
    %< >
    if i==5 %|| i==7% || i==8%&& i<=36 || i>40 && i<=44
        f(i)=l*1;
    else
       f(i)=0;
    end
   % if i==6 || i==8 % || i>36 && i<=40 || i>44 
    %  f(i)=-l*1;   
   % end
    
    %if i==41
     %   f(i)=1;
    %end    
    %if i==43
    %   f(i)=-1;
    %end    
end
%%
%Transfer function Multi forcing
a=f*V;
%x=NUM/DEN
for i=1:k
    NUM(i)=V(r,i)*a(i);%Not normalized eigenvectors, h:hysteric
end
for j=1:k
    for i=1:length(fdom)
        DEN=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        x(i,j)=NUM(j)/DEN;%Mode contribution.
    end
end
xa=abs(x);%Absolute part
xi=imag(x);%Imaginary part
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(xa(:,i)))
hold on
end
title('48 modes contribution- Multiple Forces- Hysteric Damping')
grid on
hold off

for i=1:length(fdom)
    xtc(i)=sum((x(i,:)));%Total COMPLEX transfer function.
end

for i=1:length(fdom)
    xta(i)=abs(xtc(i));%Absolute part of the total transfer function.
end
%{
for i=1:length(fdom)
    xti(i)=sum((xi(i,:)));%Phase part of the total transfer function.
end
%}
%phase
for i=1:length(fdom)
    ph(i)=angle(xtc(i));%Phase part of the total complex transfer function.
end
ph=ph*180/pi;
%%
%Different forcing plotting results

figure('units','normalized','outerposition',[0 0 1 1])

subplot(4,3,[1,2 4,5])

title('Displacement and Phase- Multiple forces- Hysteric Damping')
yyaxis left
xlabel('Frequency domain')
ylabel('Displacement')
yl=abs(xta);
semilogy(fdom,yl,'b')

yyaxis right
ylabel('Phase')
yr=ph;
plot(fdom,yr)
grid on

%{
%Checking the maximums
hold on
for j=1:length(maxp)
for i=1:length(fdom)
    if i==maxp(j)
        maxxr(i)=maxv(j);
    end
    %else
     %   maxxr(i)=0;
    %end
end
end
semilogy(fdom, maxxr,'-o')
hold off
%}
%Modes excited

%find the 2 modes that colaborate most: the ones excited.

for i=1:2*N-1
    if round(d(i)^(0.5),3)==round(d(i+1)^(0.5),3)
        xawdb(:,i)=zeros(1,length(xa(:,1)));%Xa without double modes
    else
        xawdb(:,i)=xa(:,i);
    end       
end
xawdb(:,2*N)=xa(:,2*N);

for i=1:2*N
[maxv(i), maxp(i)]=max(xawdb(:,i));%Value and position of the maximium of every mode contribution
end
[temp,originalpos] = sort( maxv, 'descend' );
n = temp(1:2);
p =originalpos(1:2);%FRF

subplot(4,3,3)
  
    stem(Vnd(:,p(1)), 'r')
    hold on
    stem(Vnb(:,p(1)))
    xlim([1 N])
    %ylim([-0.3 0.3])
    label=num2str(d(p(1))^0.5);
    title(label)
    hold off
 
 
subplot(4,3,6)
  
    stem(Vnd(:,p(2)), 'r')
    hold on
    stem(Vnb(:,p(2)))
    xlim([1 N])
    %ylim([-0.3 0.3])
    hold off 
    label=num2str(d(p(2))^0.5);
    title(label)
    legend('Disk', 'Blade')


%Force
subplot(4,3,12)
stem(f)
xlim([1 2*N])
xlabel('DOF')
title('Excitation Force')

% Modes Contribution
subplot(4,3,9)
for i=1:2*N
semilogy(fdom,abs(xa(:,i)))
hold on
end
title('48 modes contribution')
grid on
hold off

% Phase only
subplot(4,3,[7,8])
title('Displacement Multiple forces- Hysteric Damping')
plot(fdom,ph)
legend('Phase')
grid on

% Transfer Function only
subplot(4,3,[10,11])
title('Displacement Multiple forces- Hysteric Damping')
xlabel('Frequency domain')
semilogy(fdom,yl,'b')
legend('Displacement')
%ylim([ (min(yl)) (max(yl))])
grid on
%%
%Displaying values of interest
fprintf('Disk Mass : %g \n', md(1));
fprintf('Blade Mass : %g \n', mb(1));
fprintf('Response Point : %g \n',r);
fprintf('Excitation Point(for simple excitation case) : %g \n',e);
fprintf('Simple Damping Const : %g \n',b);
fprintf('Hysteric Damping Stiffness Const : %g \n',b2);
fprintf('Hysteric Damping Mass Const : %g \n',g);

