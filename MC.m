function [ msfTa, msfT, msfMa, msfM, f2T, f2M ] = MC( s,bo)
% Montecarlo study for the msf and max suppression for different multi
%forces of suppression
%% General Control Parameters:
N=4;%2*N mode shapes.
Nb=1*N;%number of blades.
%c=4;%Number of columns of mode shapes to display.
e=8;% Excitation Point
%sine=0;%Sinusoidal excitation IN THE BLADES
%sinp=1;%period of the sinusoidal exitation (ND)
%r=7;% Response Point/Drive Point
%pp=5;%Piezo position
w=1; %1 for windowing, else:0
%s=8; %Mode to be supressed with the window(up to 2*N
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
%MULT=0;%Multiple forcing yes=1, no=0;
%p=1;%percentage of the excitation force chosen to be exerted by the piezo
%bo=1;%analise blades-only in the Modal Scale Factor (1:yes, 0:no)
%% Mistuning Parameters

MISTM=1; %add Mass mistuning 0=No, 1=yes
pm=10; %percentage of change of Ma (+ve add, -ve substract)
altm=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
dofm=0; %degree of freedom to be changed on M (Nb+N Masses)
sinm=0;%sinusoisal mistuning(1 if yes, 0 if no)
sinpm=0;%period of the sinusoidal mistuning

%% Creating an excitation in any/various DOF.(Just for MIMO model)

f= zeros(1, Nb+N);
f(e)=1;
%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kb=1000*ones(N,1);%stiffness of the blade 'i'
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'
z=zeros(N,1);

dM= kdl+kdr+kb+kg;%elements of K from the equations of motion of the disks
dM= [dM;kb];%concatenate elements from equations of motion of the blades

KD= diag(dM);%create basic Stiffness matrix diagonal, without coupling elements 
Kr=diag([(-kdr(1:N-1));z],1);%to the right of the diagonal
Krc=diag((-kb), N);%to the right corner of the diagonal
KM=Kr+Krc;
KM(1,N)=-kdr(N);%term due to the left of the first disk
KM=KM+KM.';%Symmetric
KM=KM+KD;%STIFFNESS MATRIX 
KT=KM;

%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mb=1*ones(N,1);%masses of the blades
di=[md;mb];
Ma=diag(di);%MASS MATRIX

if MISTM==1
    [ MM, Ma ] = MistuningMass( Ma, pm , dofm, N, sinm, sinpm,1, altm );
else
    MM=Ma;
end
%% EIGEN PROBLEM
[VT,DT] = eig(KT,Ma);% Tuned
dT=diag(DT);
[VM,DM] = eig(KM,MM);%MisTuned
dM=diag(DM);
%% Frequency Domain
fdommacx=max(sqrt(dT(N+Nb)),sqrt(dM(N+Nb)));
fdom=0.001:0.001:1.2*(fdommacx);%Frequency domain

%% Multiple RANDOM Forces

[ weT,f2T, VminT ] = Suppressfindforce2( 0, f, VT,s );
[ weM,f2M, VminM ] = Suppressfindforce2( 0, f, VM,s );

%% Displacement of modes with and w/o mistuning

[ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g );

[ VTTT, VMMM, y2, y3 ] = suppeffect(VTT, VMM, DAMP, PHT, PHM );%here I used suppeffect function(which was created for picking up the displacement from the FRF for suppression comparison)

%% Suppression effect on the tuned model
[ VTs, VTns, PHTs, PHTns, abslambda, xrTs, xtrTs, phrTs, xrTns, xtrTns, phrTns ] = SuppComparison( 1,dT, fdom, VT, f, f2T, w, ra, s, b, g );%Damping set to 1
[ VTTs, VTTns, y2s, y3Ts ] = suppeffect(VTs, VTns, 1, PHTs, PHTns );

%% Suppression effect on the mistuned model
[ VMs, VMns, PHMs, PHMns, abslambda, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2M, w, ra, s, b, g );%Damping set to 1
[ VMMs, VMMns, y2s, y3Ms ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );

%% MSF
[ msfT ] = MSF( VTTs(:,s), VTTns(:,s), bo, Nb );
[ msfM ] = MSF( VMMs(:,s), VMMns(:,s), bo, Nb );
[ msfTa ] = MSF( VTTs(:,s), VTTns(:,s), 0, Nb );%all elements
[ msfMa ] = MSF( VMMs(:,s), VMMns(:,s), 0, Nb );%all elements

f2T=f2T';
f2M=f2M';
end