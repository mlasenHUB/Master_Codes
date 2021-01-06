function [ Vn, D ] = MSandFRF(e,r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%e: Excitation point
%r: Response point
clear all
close all

%%CONTENTS
%(1)Definition of K and M for 24 blades+disks (48 DOF)
%(2)Solution of eigen problem
%(3)Plotting mode shapes
%(4)Finding Transference function

%%
%K-STIFFNESS MATRIX
N=24;
kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
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
md=20*ones(N,1);%masses of the disks
mb=ones(N,1);%masses of the blades
d=[md;mb];
M=diag(d);%MASS MATRIX


%%
%EIGEN PROBLEM
[V,D,W] = eig(K,M);%V:eigen vectors, D:eigen values.
w=diag(D);

%%
%MASS NORMALISATION
mr=V.'*M*V;%modal mass calculation
Vn=V*mr^(0.5);%Normalised eigenvectors

for i=1:2*N
    len(i)=norm(Vn(:,i));
end
%%
%MODE SHAPES
Vnd=Vn(1:24,:);%Disks Modal Shapes
Vnb=Vn(25:2*N,:);%Blades Modal Shapes
Vndw2m=[];
Vnbw2m=[];

for i=1:2*N-1

    if round(w(i),4)~=round(w(i+1),4)
            Vndw2m=[Vndw2m Vnd(:,i)];%nomalised disk modal shapes without double modes
            Vnbw2m=[Vnbw2m Vnb(:,i)];%nomalised blade modal shapes without double modes
    end
    
end
Vndw2m=[Vndw2m  Vnd(:,48)];
Vnbw2m=[Vnbw2m  Vnb(:,48)];
%%
%PLOTTING
figure
for i=1:2*N
%i=2*a-1;
subplot(6,8,i)       % add first plot in 2 x 2 grid
stem(Vnd(:,i), 'r')
hold on
stem(Vnb(:,i))
xlim([1 24])
ylim([-0.4 0.4])
%legend('Disks', 'Blades')
hold off
end
suptitle('All Mode Shapes')

figure
for i=1:length(Vndw2m)

subplot(6,8,i)       
stem(Vndw2m(:,i), 'r')
hold on
stem(Vnbw2m(:,i))
xlim([1 24])
ylim([-0.4 0.4])
%legend('Disks', 'Blades')
hold off
end
suptitle('Without Double Modes Shapes')

%%
%TRANSFER FUNCTIONS


fdom=0:0.005:2*(max(D(:))^0.5);%Frequency domain
d=diag(D);%just the frequencies
k=length(d);

for i=1:k
    A(i)=Vn(e,i)*Vn(r,i);
end
for j=1:k
    for i=1:length(fdom)
        H(i,j)=A(j)/(d(j)-fdom(i)*fdom(i));
    end
end
H=abs(H);

%%
%PLOTTING TRANSFER FUNCTIONS

figure
for i=1:48
semilogy(fdom,H(:,i))
hold on
end
title('All 48 modes contribution')
hold off

for i=1:length(H(:,1))
    Ht(i)=sum(H(i,:));
end
figure
semilogy(fdom,Ht)
title('Total transference function')
end