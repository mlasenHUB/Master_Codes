clear all
close all

%%CONTENTS
%(1)Definition of K and M for 24 blades+disks (48 DOF)
%(2)Solution of eigen problem
%(3)Plotting of mode shapes nh
%%
%General Control Parameters:
N=4;%2*N mode shapes.
c=4;%Number of columns of mode shapes to display.
e=5;% Excitation Point
r=5;% Response Point/Drive Point
P=-2; %0.7996%0.7996(for pp=s=2);% Force applied
pp=5;%Piezo position
window=0; %1 for windowing, else:0
s=3; %Mode to be supressed with the window
ra=0.3; %Range of Piezo Forcing Frequency

%Creating an excitation in any/various DOF.(Just for MIMO model)
%Freqency domain:Total
for i=1:2*N
    %<alt 60 >alt 62
    if i==5 %&& i<=6 %|| i==2% && i<=36 || i>40 && i<=44
       f(i)=1;
    else
       f(i)=0;
    end
    if i==7 %&& i<=8 %|| i==4% && i<=40 || i>44 
       f(i)=-1;   
    end
    
    %if  i==4% && i<=40 || i>44 
    %    f(i)=-P*1;   
    %end
  %}
end
%Frequency domain: Mode to be supressed
%Piezo

for i=1:2*N
    
    if  i==pp% && i<=36 || i>40 && i<=44
        f2(i)=P*1;
    else
        f2(i)=0;
    end
end
%%
%K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kb=1000*ones(N,1);%stiffness of the disk-blade spring of disk-blade 'i'
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'
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
md=20*ones(N,1);%masses of the disks
mb=1*ones(N,1);%masses of the blades
di=[md;mb];
Ma=diag(di);%MASS MATRIX


%%
%EIGEN PROBLEM
[V,D] = eig(K,Ma);%V:eigen vectors, D:eigen values.

%THRESHOLD AT NODAL DIAMETER
for i=1:2*N
    
    maxV=max(abs(V(:,i)));
    V(:,i)=V(:,i).*(abs(V(:,i))> maxV/2);
end
%%
%MASS NORMALISATION
MR=V'*(Ma*V);%modal mass calculation
KR=V'*(K*V);%modal stiffness calculation
%Vn=V*(MR^(0.5));%Normalised eigenvectors
d=diag(D);%just the (squared) frequencies

for i=1:2*N
    len(i)=norm(V(:,i));
end
%%
%MODE SHAPES
Vnd=V(1:N,:);%Disks Modal Shapes
Vnb=V(N+1:2*N,:);%Blades Modal Shapes

%%
%PLOTTING Mode Shapes
figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(2*N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
bb=2*N/c;
for i=1:2*N

subplot(bb,c,i)       
stem(Vnd(:,i), 'r')
hold on
stem(Vnb(:,i),'-x')
xlim([1 N])
yli=-(real(max(V(:))));
yls=(real(max(V(:))));
ylim([yli yls])
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
set(get(gca,'title'),'Position',[N/2 yls 0.0])
hold off
end


figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(2*N),' Mode Shapes', ' Disks: -o ','Blades: -x','. Positive clockwise, normalise to the max. amplitude'  ])
bb=2*N/c;

for i=1:2*N

subplot(bb,c,i) 

circle(2,Vnd(:,i),Vnb(:,i));
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
hold on
end
hold off


%%
%TRANSFER FUNCTIONS


fdom=[1:0.01:2*(max(D(:))^0.5)];%Frequency domain
k=length(d);

for i=1:k
    
    A(i)=V(e,i)*V(r,i);%Modal Constant
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
% TRANSFER FUNCTIONS

for i=1:length(fdom)
    Ht(i)=sum(H(i,:));
end


%%
%(1.a)Damping Models: SIMPLE PROPORTIONAL DAMPING C=bK.

b=0.005;%proportional coefficient
kr=diag(KR);%Diagonal of modal stiffnes matrix.
mr=diag(MR);%Diagonal of modal mass matrix.
cr=b.*kr;%Diagonal of modal damping matrix.
for i=1:k
    Aspd(i)=V(e,i)*V(r,i);%Normalized eigenvectors, SPD:Simple Prop. Damp.
end
for j=1:k
    for i=1:length(fdom)
        den=complex((kr(j)-mr(j)*fdom(i)^2), fdom(i)*cr(j));%denominator.
        Hspd(i,j)=Aspd(j)/den;%Mode contribution SPD:Simple Prop. Damp.
    end
end
Hspda=abs(Hspd);%Absolute value

for i=1:length(fdom)
    Htspdc(i)=sum((Hspd(i,:)));%Total Complex transfer function.
end
for i=1:length(fdom)
    Htspda(i)=abs(Htspdc(i));%Absolute part of the total transfer function.
end

%%
%(1.b)Damping Models: PROPORTIONAL DAMPING, HYSTERIC C=H=b2K+gM
b2=0.005;%proportional coefficient- Stiffness
g=0.5;%proportional coefficient- Mass
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

for i=1:length(fdom)
    Hhtc(i)=sum((Hh(i,:)));%Complex total transfer function.
end

for i=1:length(fdom)
    Hhta(i)=abs(Hhtc(i));%Absolute part of the total transfer function.
end

%%
%(1.c)Damping Models: NON- PROPORTIONAL DAMPING ??? PENDING

%%
%%Ploting Figures with single forcing:
%(a): No Damping.
%(b): Simple Proportional Damping.
%(c): Hysteric Damping

figure('units','normalized','outerposition',[0 0 1 1])
%(a)
subplot(3,3,1)

for i=1:2*N
semilogy(fdom,abs(H(:,i)))
hold on
end
grid on
hold off
ylabel('No Damping')

subplot(3,3,2)

semilogy(fdom,abs(Ht))
grid on
xlabel('Frequency')
ylabel('x/F')

subplot(3,3,3)

plot(real(Ht), imag(Ht))
grid on
xlabel('Real x/F')
ylabel('Imaginary x/F')


%(b)

subplot(3,3,4)

for i=1:2*N
semilogy(fdom,abs(Hspda(:,i)))
hold on
end
grid on
hold off
ylabel('Simple Proportional Damping')
subplot(3,3,5)

semilogy(fdom,abs(Htspda))
grid on
xlabel('Frequency')
ylabel('x/F')

subplot(3,3,6)

plot(real(Htspdc), imag(Htspdc))
grid on
xlabel('Real x/F')
ylabel('Imaginary x/F')

%(c)

subplot(3,3,7)

for i=1:2*N
semilogy(fdom,abs(Hha(:,i)))
hold on
end
grid on
hold off
ylabel('Hysteric Damping')
subplot(3,3,8)

semilogy(fdom,abs(Hhta))
grid on
xlabel('Frequency')
ylabel('x/F')

subplot(3,3,9)

plot(real(Hhtc), imag(Hhtc))
grid on
xlabel('Real x/F')
ylabel('Imaginary x/F')

suptitle(['Excitation Point: ', num2str(e), ' Response Point: ', num2str(r), 'Single Exitation force'])

%%
%(2.a)Forcing Models: MULTIPLE FORCES, SAME MAGNITUDE. and PHASE.(Hysteric
%Damping)

%Transfer function Multi forcing
a=f*V;
%x=NUM/DEN
af=a+f2*V;%Piezo forcing in a certain frequency
for i=1:k
    
    NUM(i)=V(r,i)*a(i);%Normalized eigenvectors V
end
%Modal Constant with forcing
for i=1:k
    
    NUMf(i)=V(r,i)*af(i);%Normalized eigenvectors V
end
for j=1:k
    for i=1:length(fdom)
                
        DEN=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        x(i,j)=NUM(j)/DEN;%Mode contribution.
        %Windowing
        if window==1 && fdom(i)>=(d(s)^0.5-ra) && fdom(i)<=(d(s)^0.5+ra)% && j==s
            
           x(i,j)=NUMf(j)/DEN;%Mode contribution.
        end
        if window==0 
            
           x(i,j)=NUMf(j)/DEN;%Mode contribution.
        end
    end
end

xa=abs(x);%Absolute part
xi=imag(x);%Imaginary part
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:2*N
semilogy(fdom,abs(xa(:,i)))
[maxv(i), maxp(i)]=max(xa(:,i));%Value and position of the maximium of every mode contribution
hold on
end

title(['Number of DOF(Modes): ', num2str(2*N), ' -Multiple Forces (Excitation and Piezo)-; Piezo pos: ', num2str(pp), ' - Hysteric Damping'])
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

title('Operation Deflexion Shape and Phase- Multiple forces- Hysteric Damping')
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

%Modes excited

%find the 2 modes that colaborate most: the ones excited.


for i=1:2*N
[maxv(i), maxp(i)]=max(xa(:,i));%Value and position of the maximium of every mode contribution
end
[temp,originalpos] = sort( maxv, 'descend' );
n = temp(1:2);
p =originalpos(1:2);%FRF

subplot(4,3,3)
  
    stem(Vnd(:,p(1)), 'r')
    hold on
    stem(Vnb(:,p(1)))
    xlim([1 N])
    ylim([-0.4 0.4])
    label=num2str(d(p(1))^0.5);
    title(label)
    hold off
 
 
subplot(4,3,6)
  
    stem(Vnd(:,p(2)), 'r')
    hold on
    stem(Vnb(:,p(2)))
    xlim([1 N])
    ylim([-0.4 0.4])
    hold off 
    label=num2str(d(p(2))^0.5);
    title(label)


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
for i=1:2*N
semilogy(fdom,abs(xa(:,i)))
hold on
end
title([num2str(2*N), ' modes contribution'])
grid on
hold off

% Phase only
subplot(4,3,[7,8])
plot(fdom,ph)
legend('Phase')
grid on


% Transfer Function only
subplot(4,3,[10,11])
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

%%
%Finding the best positions, forcing values to suppress certain modes

%F=1;
%Mmin

for k = 1:2*N
    f2=zeros(1, 2*N);
    for ii =1:2*N
        f2(ii)=1;
        a2=f2*V;
        M(ii,k)= a(k)+a2(ii);%k: mode, ii:force position.(a(k):'pseudo' modal constant
    end   
end

%Best positions to place the piezo to suppress that mode:
Mmin=zeros(4, 2*N);
for i =1:2*N
    [minval, pos]=min(abs(M(:,i)));%We want to minimize the final 'Pseudo Modal constant'
    Mmin(1, i)= pos;
    Mmin(2, i)= M(pos, i);
    Mmin(3, i)= min(a(i));
    if abs(Mmin(2,i)) < abs(Mmin(3,i))
        Mmin(4, i)= 1;%1 is I reduced it, else 0
    end
    Mmin(5, i)=a2(pos);
   
end
Mnom=M./(max(abs(M(:))));

figure('units','normalized','outerposition',[0 0 1 1])

for k=1:2*N
    
subplot(bb+1,c,k)

plot(Mnom(:,k))
hold on
plot(abs(Mnom(:,k)), 'r')
hold on
plot(a(k),'--gx')
hold on
plot(a2,'--ko')
hold off
end
legend('Normalise Pseudo Modal Constant(PMC)', 'Positive PMC','Exitation contribution to PMC','Effect of piezo position in PMC')
subplot(bb+1,c,[2*N+1 2*N+c])
for i=1:2*N
    plot(Mnom(:,i))
    hold on
end
hold off
%Mmin2
%Best forcing value to place it at any position, to suppress that mode
Mmin2=zeros(4, 2*N);
for k = 1:2*N%Mode Shape
    for ii =1:2*N%Forcing pos
        Mmin2(ii, k)= -a(k)/V(ii, k);
    end    
end
Mmin2(~isfinite(Mmin2))=NaN;%Neglect infinite values of forcing
Mmin2(isnan(Mmin2))=0;
Mnom2=Mmin2./max(abs(Mmin2(:)));
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

%% DISPLAY DISPLACEMENT OF MODES WITH PIEZO SUPPRESSION w7o Damping
figure('units','normalized','outerposition',[0 0 1 1])
suptitle(['Mode Supressed: ', num2str(s),'No-Dumping', ' Disks: -o ','Blades: -x' ])
for i=1:2*N%Go through all the modes to see the effect of the piezo

[val indf]=min(abs((d(i)^0.5-fdom)));%Find the index of the nat.freq. of the mode to suppress in the freq. domain.

Vsup=H(indf, :);
Vs(:,i)=Vsup;
subplot(bb,c,i) 
circle(2,Vsup(1:N),Vsup(N+1:2*N));
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
hold on
end
hold off

%% SUPPRESSION FACTOR
%Local Supression Factor: SFl
for i=1:2*N
    SFl(i)=abs(norm(V(:,i)-Vs(:,i)))/norm(V(:,i));
end

%Total Suppression factor: SF (Effect of the piezo over all the modes)
SFnum=0;
SFden=0;
for i=1:2*N
     
    SFnum=abs(norm(V(:,i)-Vs(:,i)))+SFnum;
    SFden=norm(V(:,i))+SFden;
 end
 SF= SFnum/SFden;
 
 figure
 plot(SFl)
 xlabel('Modes')
ylabel('Local Suppression Factor')




