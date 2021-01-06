%% Model 3: 1 Mass disk + 2 masses blades in parallel.
clear all
close all
%% General Control Parameters:
N=4;% number of disks
Nb=2*N;%Number of blades per disk
c=4;%Number of columns of mode shapes to display.
e=5;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=2;%period of the sinusoidal exitation (ND)
r=7;% Response Point/Drive Point
pp=5;%Piezo position
w=1; %1 for windowing, else:0
s=8; %Mode to be supressed with the window(up to 2*N
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=1;%1 for plotting, 0 for no plotting
t=1; %no torsional modes in 2nd model
%% Mistuning Parameters

MISTK=0; % add stiffness 0=No, 1=yes 
pk= 10;%percentage of change of K (+ve add, -ve substract)
dofk= 1; %degree of freedom to be changed on K (N up, down blades)
upbk=1;%1 if the mistuning is applied to the UPper mass blade (0 for the DOWN)
sinmk=0;%sinusoisal mistuning
sinpk=0;%period of the sinusoidal mistuning

MISTM=1; %add Mass mistuning 0=No, 1=yes
pm=70; %percentage of change of Ma (+ve add, -ve substract)
dofm=8; %degree of freedom to be changed on M (Nb+N Masses)
upbm=1;% 1 if the mistuning is applied to the UPper mass blade (0 for the DOWN)
sinm=1;%sinusoisal mistuning
sinpm=1;%period of the sinusoidal mistuning

%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kbu=1000*ones(N,1);%stiffness of the blade 'i' UP
kbd=1000*ones(N,1);%stiffness of the blade 'i' DOWN
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'

if MISTK==1
    if upbk==1
    [ kbuM, kbu ] = MistuningStiffness( kbu, pk, dofk, N, sinmk, sinpk  );
    kbdM=kbd;
    else
    [ kbdM, kbd ] = MistuningStiffness( kbd, pk, dofk, N, sinmk, sinpk );
    kbuM=kbu;
    end 
else
    kbuM=kbu;
    kbdM=kbd;
end
%Tuned
H= diag(kdl+kdr+kbu+kbd+kg);%elements of K from the equations of motion of the disks
H2= diag(-kdr(1:N-1),1);%concatenate elements from equations of motion of the blades
H=H+H2;
H(1,N)=-kdr(N);
H=H+triu(H,1)';

D0=diag(zeros(1,N));
D1=diag(-kbu);
D2=diag(-kbd);

KT= [H D1 D2; D1 -D1 D0; D2 D0 -D2];
%Mistuned
HM= diag(kdl+kdr+kbuM+kbdM+kg);%elements of K from the equations of motion of the disks
H2M= diag(-kdr(1:N-1),1);%concatenate elements from equations of motion of the blades
HM=HM+H2M;
HM(1,N)=-kdr(N);
HM=HM+triu(HM,1)';

D0=diag(zeros(1,N));
D1M=diag(-kbuM);
D2M=diag(-kbdM);

KM= [HM D1M D2M; D1M -D1M D0; D2M D0 -D2M];


%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mbu=2*ones(N,1);%masses of the blades LEFT
mbd=1*ones(N,1);%masses of the blades RIGHT
di=[md;mbu;mbd];
MT=diag(di);%MASS MATRIX

if MISTM==1
    [ MM, MT ] = MistuningMass( MT, pm , dofm, N, sinm, sinpm,upbm );
else 
    MM=MT;
end 
%% EIGEN PROBLEM
[VT,DT] = eig(KT,MT);% Tuned
dT=diag(DT);
[VM,DM] = eig(KM,MM);%MisTuned
dM=diag(DM);
%% Creating an excitation in any/various DOF.(Just for MIMO model)
%Freqency domain:Total
for i=1:Nb+N
    %<alt 60 >alt 62
    if i==e %|| i==e+2% && i<=36 || i>40 && i<=44
       f(i)=1;
    
    else
       f(i)=0;
    end
end

%% Finding the best positions, forcing values to suppress certain modes

%First Approach : Mmin: find how much I can suppress a mode with a given force of magnitude
F=1;
[Mmin]=Suppressgivenforce( c,f, VT,F );

%Second Approach : Mmin2: What is the best force magnitude, to place at any position, to suppress a certain mode
%Tuned
[Mmin2T, Mnom2T] = Suppressfindforce( f, VT );

Mmin2Tc=Mmin2T.*(abs(Mmin2T)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)
PT= Mmin2T(pp,s);%See matrix Mmin2% Force applied
f2T=zeros(1, Nb+N);
f2T(pp)=PT;

%Mistuned
[Mmin2M, Mnom2M] = Suppressfindforce( f, VM );
Mmin2Mc=Mmin2M.*(abs(Mmin2M)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)
PM= Mmin2M(pp,s);%See matrix Mmin2% Force applied
f2M=zeros(1, Nb+N);
f2M(pp)=PM;

%% MODE SHAPES
VdM=VM(1:N,:);%Disks Modal Shapes with mistuning
VbM=VM(N+1:end,:);%Blades Modal Shapes with mistuning

%Modes without mistuning
VdT=VT(1:N,:);%Disks Modal Shapes without Mistuning(Tuned)
VbT=VT(N+1:end,:);


%% PLOTTING Mode Shapes and Natural Frequencies vs Nodal Diameters
if plotmodes==1
Plotmodes(c, dM, VM,VdM, VbM, t )
Plotmodes(c, dT, VT,VdT, VbT, t )
end
%% Frequency Domain 
fdom=0.001:0.001:2*(max(DT(:))^0.5);%Frequency domain

%% Displacement of modes with and w/o mistuning


[ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g );
PlotTFpoints( c, xrT, xtrT,phrT, fdom, 'TUNED')  
PlotTFpoints( c, xrM, xtrM, phrM, fdom, 'MISTUNED')  

[ VTTT, VMMM, y2, y3 ] = suppeffect(VTT, VMM, DAMP, PHT, PHM );%here I used suppeffect function(which was created for picking up the displacement from the FRF for suppression comparison)
barcomparison( c, VTTT, VMMM, 'Tuned', 'Mistuned');
%% Suppression effect on the mistuned model
[ VMs, VMns, PHMs, PHMns, abslambda, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2M, w, ra, s, b, g );%Damping set to 1
PlotTFpoints( c, xrMs, xtrMs, phrMs, fdom, 'MISTUNED-SUPPRESSED')
PlotTFpoints( c, xrMns, xtrMns, phrMns, fdom, 'MISTUNED-UNSUPPRESSED')

[ VMMs, VMMns, y2s, y3s ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );
barcomparison( c, VMMs, VMMns, 'Mistuned-Suppressed', 'Mistuned-Unsuppressed');

%% General results: ratio of diplacement, MACs, natural frequencies and suppression

[MACe, MACf, r, nfT, nfM ] = GeneralandMAC( N,Nb, f, f2M, VM, VT, VMMM, VTTT, KM, KT, MM, MT, dT, dM, VMMs, VMMns );