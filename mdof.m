clear all
close all

%%CONTENTS
%(1)Definition of K and M for 24 blades+disks (48 DOF)
%(2)Solution of eigen problem
%(3)Plotting of mode shapes nh

%%
%K-STIFFNESS MATRIX
N=24;
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
md=10*ones(N,1);%masses of the disks
mb=10*ones(N,1);%masses of the blades
di=[md;mb];
M=diag(di);%MASS MATRIX


%%
%EIGEN PROBLEM
[V,D] = eig(K,M);%V:eigen vectors, D:eigen values.

%%
%MASS NORMALISATION
MR=V'*(M*V);%modal mass calculation
KR=V'*(K*V);%modal stiffness calculation
Vn=V*MR^(0.5);%Normalised eigenvectors
d=diag(D);%just the frequencies

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
subplot(6,8,i)       % add first plot in 2 x 2 grid
stem(Vnd(:,i), 'r')
hold on
stem(Vnb(:,i))
xlim([1 N])
ylim([-0.3 0.3])
%legend('Disks', 'Blades')
label=num2str(d(i)^0.5);
title(label)
set(get(gca,'title'),'Position',[5.5 0.2 0.5])
hold off
end
suptitle('All Mode Shapes')

%%
%TRANSFER FUNCTIONS

e=1;% Excitation Point
r=25;% Response Point
fdom=[6:0.01:2*(max(D(:))^0.5)];%Frequency domain
k=length(d);

for i=1:k
    A(i)=Vn(e,i)*Vn(r,i);%Modal Constant(V used here because they seem to be normalised)
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
Hspdr=real(Hspd);%Real part

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(Hspdr(:,i)))
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
    Htspdr(i)=sum((Hspdr(i,:)));%Real part of the total transfer function.
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
semilogy(fdom,abs(Htspdr))
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
Hhr=real(Hh);%Real part

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(Hhr(:,i)))
hold on
end
title('All 48 modes contribution with Hysteric Damping')
grid on
hold off
for i=1:length(fdom)
    Hhtc(i)=sum((Hh(i,:)));%Complex total transfer function.
end

for i=1:length(fdom)
    Hhtr(i)=sum((Hhr(i,:)));%Real part of the total transfer function.
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
semilogy(fdom,abs(Hhtr))
grid on

hold off
subplot(1,2,2)
plot(real(Hhtc), imag(Hhtc))
grid on
suptitle('Total transference function Hysteric Damping')
%%
%(1.c)Damping Models: NON- PROPORTIONAL DAMPING ???

%%
%(2.a)Forcing Models: MULTIPLE FORCES, SAME MAGNITUDE. and PHASE.(Hysteric
%Damping)
%Creating an excitation in any/various DOF.
for i=1:2*N
    %< >
    if i>=25 && i<=28 || i>32 && i<=36 || i>40 && i<=44
        f(i)=1;
     else
        f(i)=1;
    end
    if i>28 && i<=32 || i>36 && i<=40 || i>44 
        f(i)=1;   
    end
    %if i==41
     %   f(i)=1;
    %end    
    %if i==43
    %   f(i)=-1;
    %end    
end
%%
%Transfer function Multi forcing
a=f*Vn;
%x=NUM/DEN
for i=1:k
    NUM(i)=Vn(r,i)*a(i);%Not normalized eigenvectors, h:hysteric
end
for j=1:k
    for i=1:length(fdom)
        DEN=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        x(i,j)=NUM(j)/DEN;%Mode contribution.
    end
end
xr=real(x);%Real part
xi=imag(x);%Imaginary part
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(xr(:,i)))
[maxv(i), maxp(i)]=max(xr(:,i));%Value and position of the maximium of every mode contribution
hold on
end
title('48 modes contribution- Multiple Forces- Hysteric Damping')
grid on
hold off

for i=1:length(fdom)
    xtr(i)=sum((xr(i,:)));%Real part of the total transfer function.
end

for i=1:length(fdom)
    xti(i)=sum((xi(i,:)));%Phase part of the total transfer function.
end
%phase
for i=1:length(fdom)
    ph(i)=atan(abs(xti(i)/xtr(i)));%Phase part of the total transfer function.
end
%%
%Different forcing plotting results: Force 1
%find the 2 modes that colaborate most: the ones excited.
[temp,originalpos] = sort( maxv, 'descend' );
n = temp(1:2);
p =originalpos(1:2);%FRF
figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,3,[1,6])
semilogy(fdom,abs(xtr))
grid on
title('Displacement- Multiple forces- Hysteric Damping')
xlabel('Frequency domain')
ylabel('Displacement')

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
subplot(4,3,7)
  
    stem(Vnd(:,p(1)), 'r')
    hold on
    stem(Vnb(:,p(1)))
    xlim([1 N])
    ylim([-0.3 0.3])
    label=num2str(d(p(1))^0.5);
    title(label)
    hold off
    
    
subplot(4,3,10)
  
    stem(Vnd(:,p(2)), 'r')
    hold on
    stem(Vnb(:,p(2)))
    xlim([1 N])
    ylim([-0.3 0.3])
    hold off 
    label=num2str(d(p(2))^0.5);
    title(label)


%Force
subplot(4,3,[8,9])
stem(f)
xlim([1 2*N])
xlabel('DOF')
title('Excitation Force')

%$8 Modes Contribution
subplot(4,3,[11,12])
for i=1:2*N
semilogy(fdom,abs(xr(:,i)))
hold on
end
title('48 modes contribution- Multiple Forces- Hysteric Damping')
grid on
hold off
%%
%Displaying values of interest
fprintf('Disk Mass : %g \n', md(1));
fprintf('Blade Mass : %g \n', mb(1));
fprintf('Response Point : %g \n',r);
fprintf('Excitation Point(for simple excitation case) : %g \n',e);
fprintf('Simple Damping Const : %g \n',b);
fprintf('Hysteric Damping Stiffness Const : %g \n',b2);
fprintf('Hysteric Damping Mass Const : %g \n',g);

