%Study case
clear all
close all
%% General Control Parameters:
N=4;% number of disks
Nb=2*N;%Number of blades per disk
c=4;%Number of columns of mode shapes to display.
e=12;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=2;%period of the sinusoidal exitation (ND)
%r=7;% Response Point/Drive Point
DOUBLE=1; % DOUBLE forcing activated(1: yes, 0:no)
ppm=2;% piezos per mode to suppressed (1=1 piezo per mode, 2=2 piezos per mode)
w=1; %1 for windowing, else:0
s=1; %Mode to be supressed with the window(up to Nb+N)
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=1;%1 for plotting, 0 for no plotting
t=0; %no torsional modes in 2nd model
p=1;%percentage of the excitation force chosen to be exerted by the piezo
%% Mistuning Parameters

MISTK=1; % add stiffness 0=No, 1=yes 
pk= 10;%percentage of change of K (+ve add, -ve substract)
dofk= 1; %degree of freedom to be changed on K (N left, right blades)
leftbk=1;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinmk=1;%sinusoisal mistuning
sinpk=1;%period of the sinusoidal mistuning

MISTM=0; %add Mass mistuning 0=No, 1=yes
pm=70; %percentage of change of Ma (+ve add, -ve substract)
dofm=8; %degree of freedom to be changed on M (Nb+N Masses)
leftbm=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinm=1;%sinusoisal mistuning(1 if yes, 0 if no)
sinpm=1;%period of the sinusoidal mistuning

%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kbl=1000*ones(N,1);%stiffness of the blade 'i' LEFT
kbr=1000*ones(N,1);%stiffness of the blade 'i' RIGHT
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'

if MISTK==1
    if leftbk==1
    [ kblM, kbl ] = MistuningStiffness( kbl, pk, dofk, N, sinmk, sinpk  );
    kbrM=kbr;
    else
    [ kbrM, kbr ] = MistuningStiffness( kbr, pk, dofk, N, sinm, sinpm  );
    kblM=kbl;
    end 
else
    kblM=kbl;
    kbrM=kbr;
end
%Tuned
H= diag(kdl+kdr+kbl+kg);%elements of K from the equations of motion of the disks
H2= diag(-kdr(1:N-1),1);%concatenate elements from equations of motion of the blades
H=H+H2;
H(1,N)=-kdr(N);
H=H+triu(H,1)';

D0=diag(zeros(1,N));
D1=diag(-kbl);
D2=diag(kbl+kbr);
D3=diag(-kbr);
D4=diag(kbr);

KT= [H D1 D0; D1 D2 D3; D0 D3 D4];
%Mistuned
HM= diag(kdl+kdr+kblM+kg);%elements of K from the equations of motion of the disks
H2M= diag(-kdr(1:N-1),1);%concatenate elements from equations of motion of the blades
HM=HM+H2M;
HM(1,N)=-kdr(N);
HM=HM+triu(HM,1)';

D1M=diag(-kblM);
D2M=diag(kblM+kbrM);
D3M=diag(-kbrM);
D4M=diag(kbrM);

KM= [HM D1M D0; D1M D2M D3M; D0 D3M D4M];


%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mbl=2*ones(N,1);%masses of the blades LEFT
mbr=1*ones(N,1);%masses of the blades RIGHT
di=[md;mbl;mbr];
Ma=diag(di);%MASS MATRIX

if MISTM==1
    [ MM, Maor ] = MistuningMass( Ma, pm , dofm, N, sinm, sinpm,leftbm );
else 
    MM=Ma;
end 
%% EIGEN PROBLEM
[VT,DT] = eig(KT,Ma);% Tuned
dT=diag(DT);
[VM,DM] = eig(KM,MM);%MisTuned
dM=diag(DM);
%% Creating an excitation in any/various DOF.(Just for MIMO model)
%Freqency domain:Total
f= zeros(1, Nb+N);
f(e)=1;

%% Finding the best positions, forcing values to suppress certain modes

%TYPE A: GIVEN A FORCE , HOW MUCH I CAN SUPPRESS
nrT=b+g./dT;
nrM=b+g./dM;
%Approach 3: 2 PIEZOS/2DOF
%Tuned
[ MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T ] = SuppressgivenforceA3( VT,nrT,p,dT,f,N,Nb);
%Mistuned
[ MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M ] = SuppressgivenforceA3( VM,nrM,p,dM,f,N,Nb);

%% Optimization of the LOCATION of the suppression force for minimum MSF(FORCE is CHOSEN)
%Type A optimization: Minimum MSF

if DOUBLE==1 && ppm==2
    [ PPT, MSFfT, idlocT1A3, idlocT2A3, OPTT ] = OptPiezoA3(MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T,s, N,Nb);
    [ PPM, MSFfM, idlocM1A3, idlocM2A3, OPTM ] = OptPiezoA3(MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M,s, N,Nb);
    
end

%% Creation of the Optimum suppression force (TYPE B)
if DOUBLE==1 && ppm==2%TWO piezo per mode
    [ f2M ] = CreateSuppForce(MminB3T, MminB3M , MISTK, MISTM,idloc, Pt, N,Nb,s,ppm);%create suppressing force.
end

%% Frequency Domain 
fdom=0.001:0.001:1.2*(max(DT(:))^0.5);%Frequency domain