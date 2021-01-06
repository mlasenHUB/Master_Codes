%% Model 2: 1 Mass disk + 2 masses blades. Just Tuned and Mistuned Torsion (No linear)
clear all
close all
%% General Control Parameters:
N=4;% number of disks
Nb=2*N;%Number of blades per disk
c=4;%Number of columns of mode shapes to display.
e=12;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=2;%period of the sinusoidal exitation (ND)
r=7;% Response Point/Drive Point
pp=5;%Piezo position(blisk unit)
DOUBLE=1; % DOUBLE forcing activated(1: yes, 0:no)
ppd=1;%Piezo postition DOUBLE(0: disk+lb, 1: lb+rb, 2:rb+disk)
w=1; %1 for windowing, else:0
s=8; %Mode to be supressed with the window(up to Nb+N)
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=1;%1 for plotting, 0 for no plotting
t=0; %no torsional modes in 2nd model
%% Mistuning Parameters

MISTK=0; % add stiffness 0=No, 1=yes 
pk= 10;%percentage of change of K (+ve add, -ve substract)
dofk= 1; %degree of freedom to be changed on K (N left, right blades)
leftbk=1;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinmk=0;%sinusoisal mistuning
sinpk=0;%period of the sinusoidal mistuning

MISTM=1; %add Mass mistuning 0=No, 1=yes
pm=70; %percentage of change of Ma (+ve add, -ve substract)
dofm=8; %degree of freedom to be changed on M (Nb+N Masses)
leftbm=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinm=1;%sinusoisal mistuning(1 if yes, 0 if no)
sinpm=1;%period of the sinusoidal mistuning

%% K-STIFFNESS MATRIX

kdr=1000000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kbl=100*ones(N,1);%stiffness of the blade 'i' LEFT
kbr=100*ones(N,1);%stiffness of the blade 'i' RIGHT
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


%% M- Inertia-Matrix
md=30*ones(N,1);%masses of the disks
mbl=2*ones(N,1);%masses of the blades LEFT
mbr=1*ones(N,1);%masses of the blades RIGHT
%mass to inertia
rad=2;%radious of the plate representing the blades.
md=md*((1000*rad)^2/4);%masses of the disks
mbl=mbl*(rad^2/4);%inertia of the LEFT plate
mbr=mbr*(rad^2/4);%inertia of the RIGHT plate

di=[md;mbl;mbr];
MT=diag(di);%MASS MATRIX

if MISTM==1
    [ MM, MT ] = MistuningMass( MT, pm , dofm, N, sinm, sinpm,leftbm );
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

%% Finding the best positions, forcing values to suppress certain modes (various approaches)

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


%Third Approach: suppressing force acting in 2DOF per blisk unit
%Tuned
[Mmin3T, Mnom3T] = Suppressfinddoubleforce( f, VT, 2 );
%Mistuned
[Mmin3M, Mnom3M] = Suppressfinddoubleforce( f, VM, 2 );

if DOUBLE==1
    Pd= Mmin3T(pp,s);
    if MISTK ==1 || MISTM==1 %if there is mistuning use the mistuned optimised matrix
        Pd = Mmin3M(pp, s);
    end
    f2d=zeros(1, Nb+N);
    f2d(pp)=Pd;
    
    if ppd==0 || ppd==1
        f2d(pp+N)=-Pd;
    else
        f2d(pp-2*N)=-Pd;
    end
    
f2M=f2d;
end

%Fourth approach: 2 piezos per mode
%Tuned
[Mmin4T, Mnom4T] = Suppressfinddoubleforce2( f, VT, Nb/N );

%Mistuned
[Mmin4M, Mnom4M] = Suppressfinddoubleforce2( f, VM, Nb/N );

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

%% Optimization of the amount of piezos
%1 piezo to suppress a mode
[ FtT, PtT ] = OptPiezo( Mmin3T );
[ FtM, PtM ] = OptPiezo( Mmin3M );
%2 piezos to suppress a mode
[ MmT, MmT2, PtT2,PT ] = OptPiezo2( Mmin4T );
[ MmM, MmM2, PtM2,PM ] = OptPiezo2( Mmin4M );
