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
P=0; %0.7996%0.7996(for pp=s=2);% Force applied
pp=2;%Piezo position
window=0; %1 for windowing, else:0
s=2; %Mode to be supressed with the window
ra=0.1; %Range of Piezo Forcing Frequency

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
[Vx,Dx] = eig(K,Ma);
%THRESHOLD AT NODAL DIAMETER PER EACH MODE
for i=1:2*N
    
    maxV=max(abs(Vx(:,i)));
    Vxa(:,i)=Vx(:,i).*(abs(Vx(:,i))>maxV/2);
end
mp.Digits(34);
K=mp(K);
Ma=mp(Ma);
[V,D] = eig(K,Ma,'qr');%V:eigen vectors, D:eigen values.(qr, dc, mr)
%[V1,D] = eig(K,M, 'chol');


%%
%MASS NORMALISATION
%MR=V'*(Ma*V);%modal mass calculation
%KR=V'*(K*V);%modal stiffness calculation
%Vn=V*(MR^(0.5));%Normalised eigenvectors
d=diag(D);%just the (squared) frequencies
dx=diag(Dx);%just the (squared) frequencies

for i=1:2*N
    len(i)=norm(V(:,i));
end
%%
%MODE SHAPES
Vnd=V(1:N,:);%Disks Modal Shapes
Vnb=V(N+1:2*N,:);%Blades Modal Shapes
Vndx=Vx(1:N,:);%Disks Modal Shapes w/o precision
Vnbx=Vx(N+1:2*N,:);%Blades Modal Shapes w/0 precision
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
%yli=-(real(max(V(:))));
%yls=(real(max(V(:))));
%ylim([yli yls])
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
%set(get(gca,'title'),'Position',[N/2 yls 0.0])
hold off
end

figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(2*N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
bb=2*N/c;
for i=1:2*N

subplot(bb,c,i)       
stem(Vndx(:,i), 'r')
hold on
stem(Vnbx(:,i),'-x')
xlim([1 N])
%yli=-(real(max(V(:))));
%yls=(real(max(V(:))));
%ylim([yli yls])
label=(['F', num2str(dx(i)^0.5), ' P', num2str(i)]);
title(label)
%set(get(gca,'title'),'Position',[N/2 yls 0.0])
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

figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(2*N),' Mode Shapes', ' Disks: -o ','Blades: -x','. Positive clockwise, normalise to the max. amplitude'  ])
bb=2*N/c;

for i=1:2*N

subplot(bb,c,i)       
circle(2,Vndx(:,i),Vnbx(:,i));
label=(['F', num2str(dx(i)^0.5), ' P', num2str(i)]);
title(label)
hold on
end
hold off