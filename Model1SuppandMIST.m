clear all
close all
%% General Control Parameters:
N=4;%2*N mode shapes.
Nb=1*N;%number of blades.
c=4;%Number of columns of mode shapes to display.
e=8;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=1;%period of the sinusoidal exitation (ND)
r=7;% Response Point/Drive Point
pp=5;%Piezo position
w=1; %1 for windowing, else:0
s=8; %Mode to be supressed with the window(up to 2*N
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=0;%1 for plotting, 0 for no plotting
plotdeandmc=1;%1 for plotting damping effect and modal contribution(0 for no)
t=0; %no torsional modes in 2nd model
p=1;%percentage of the excitation force chosen to be exerted by the piezo
loadopt=1; %(1 if yes)load the optimum from a previously saved mode
td= 0; %tip to disk piezo: 1=operating, 0=not in operation
bo=1;%analise blades-only in the Modal Scale Factor (1:yes, 0:no)
%% Mistuning Parameters

MISTK=0; % add stiffness 0=No, 1=yes 
pk= 10;%percentage of change of K (+ve add, -ve substract)
altk=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
dofk= nan; %degree of freedom to be changed on K (N left, right blades)
sinmk=0;%sinusoisal mistuning
sinpk=0;%period of the sinusoidal mistuning

MISTM=1; %add Mass mistuning 0=No, 1=yes
pm=10; %percentage of change of Ma (+ve add, -ve substract)
altm=0;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
dofm=5; %degree of freedom to be changed on M (Nb+N Masses)
sinm=0;%sinusoisal mistuning(1 if yes, 0 if no)
sinpm=0;%period of the sinusoidal mistuning

%% Creating an excitation in any/various DOF.(Just for MIMO model)

f= zeros(1, Nb+N);
f(e)=12;
%f(9)=1;
%f(11)=1;
%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kb=1000*ones(N,1);%stiffness of the blade 'i'
%Add mistuning in a blade
if MISTK==1
    [ kb, kbor ] = MistuningStiffness( kb, pk, dofk, N, sinmk, sinpk, altk  );
else
    kbor=kb;
end

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



%K without mistuning

dM= kdl+kdr+kbor+kg;%elements of K from the equations of motion of the disks
dM= [dM;kbor];%concatenate elements from equations of motion of the blades

KD= diag(dM);%create basic Stiffness matrix diagonal, without coupling elements 
Kr=diag([(-kdr(1:N-1));z],1);%to the right of the diagonal
Krc=diag((-kbor), N);%to the right corner of the diagonal
KT=Kr+Krc;
KT(1,N)=-kdr(N);%term due to the left of the first disk
KT=KT+KT.';%Symmetric
KT=KT+KD;%STIFFNESS MATRIX 


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

%% Finding the best positions, forcing values to suppress certain modes

%TYPE A: GIVEN A FORCE , HOW MUCH I CAN SUPPRESS
% Approach 1 : 1 PIEZO/1 DOF

%Tuned
nrT=b+g./dT;
[ MA1T, MpA1T, MnA1T ] = Suppressgivenforce( VT,nrT,p,dT,f);
%Mistuned
nrM=b+g./dM;
[ MA1M, MpA1M, MnA1M ] = Suppressgivenforce( VM,nrM,p,dM,f);

%TYPE B: FIND THE MINIMUM FORCE TO SUPPRESS A MODE COMPLETELY
%recover the displacement at resonance from the FRF.


[ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g );

if DAMP==1 %If there is DAMPING then consider the phase to see the direction of the displacement
  VTT=VTT.*sign(sin(deg2rad(PHT)));
  VMM=VMM.*sign(sin(deg2rad(PHM)));
end


% Approach 1 : 1 PIEZO/1 DOF
%Tuned
[MminB1T, MnomB1T] = Suppressfindforce( f, VTT );
MminB1Tc=MminB1T.*(abs(MminB1T)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)

%Mistuned
[MminB1M, MnomB1M] = Suppressfindforce( f, VMM );
MminB1Mc=MminB1M.*(abs(MminB1M)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)

%% Creation of the Optimum suppression force TYPE B(from the minimum)

%Just comparing the suppressing effect on the model with mistuning (neglect
%the tuned model)

%Tuned
f2T=zeros(1,N+Nb);
auxT=abs(MminB1T(:,s));
auxT(e)=nan;%do not consider force opposite in phase.
[mvT, idlocT]=min(auxT);
f2T(idlocT)=MminB1T(idlocT,s);%create suppressing force.

%Mistuned
f2M=zeros(1,N+Nb);
auxM=abs(MminB1M(:,s));
auxM(e)=nan;%do not consider force opposite in phase.
[mvM, idlocM]=min(auxM);
f2M(idlocM)=MminB1M(idlocM,s);%create suppressing force

%% MODE SHAPES
VdM=VM(1:N,:);%Disks Modal Shapes with mistuning
VbM=VM(N+1:2*N,:);%Blades Modal Shapes with mistuning

%Modes without mistuning
VdT=VT(1:N,:);%Disks Modal Shapes without Mistuning(Tuned)
VbT=VT(N+1:2*N,:);

%% PLOTTING Mode Shapes and Natural Frequencies vs Nodal Diameters
if plotmodes==1
Plotmodes(c, dM, VM,VdM, VbM, t )
Plotmodes(c, dT, VT,VdT, VbT, t )
end


%% PLOTTING DAMPING EFFECT AND MODES CONTRIBUTION AT 1 LOCATION
if plotdeandmc==1
    
    %tuned
    [ x, xtu, phu ] = NoDamp( r, VT,dT,  f, f2T, w, 0,s, fdom );%undamped
    [ xd, xtd, phd ] = ModalDamp( b,g, r, VT,dT,  f, f2T,w, 0, s, fdom );%damped
    
    PlotTF( xd,xtd, f2T, e, fdom, r )%Plot transference function and its modes contributions.
    
    PlotMain(x,xtu,xtd, phu,phd,dT, fdom, VdT, VbT, f, f2T)%Plot the summary of results: Total transference function, peaks
    
end

%% Displacement of modes with and w/o mistuning


[ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g );
PlotTFpoints( c, xrT, xtrT,phrT, fdom, 'TUNED')  
PlotTFpoints( c, xrM, xtrM, phrM, fdom, 'MISTUNED')  

[ VTTT, VMMM, y2, y3 ] = suppeffect(VTT, VMM, DAMP, PHT, PHM );%here I used suppeffect function(which was created for picking up the displacement from the FRF for suppression comparison)
barcomparison( c, VTTT, VMMM, 'Tuned', 'Mistuned');

%% Suppression effect on the tuned model
[ VTs, VTns, PHTs, PHTns, abslambdaT, xrTs, xtrTs, phrTs, xrTns, xtrTns, phrTns ] = SuppComparison( 1,dT, fdom, VT, f, f2T, w, ra, s, b, g );%Damping set to 1
PlotTFpoints( c, xrTs, xtrTs, phrTs, fdom, 'TUNED-SUPPRESSED')
PlotTFpoints( c, xrTns, xtrTns, phrTns, fdom, 'TUNED-UNSUPPRESSED')

[ VTTs, VTTns, y2Ts, y3Ts ] = suppeffect(VTs, VTns, 1, PHTs, PHTns );
barcomparison( c, VTTs, VTTns, 'Tuned-Suppressed', 'Tuned-Unsuppressed');




%% Suppression effect on the mistuned model
[ VMs, VMns, PHMs, PHMns, abslambdaM, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2M, w, ra, s, b, g );%Damping set to 1
PlotTFpoints( c, xrMs, xtrMs, phrMs, fdom, 'MISTUNED-SUPPRESSED')
PlotTFpoints( c, xrMns, xtrMns, phrMns, fdom, 'MISTUNED-UNSUPPRESSED')

[ VMMs, VMMns, y2Ms, y3Ms ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );
barcomparison( c, VMMs, VMMns, 'Mistuned-Suppressed', 'Mistuned-Unsuppressed');

%% General results: ratio of diplacement, MACs, natural frequencies and suppression

[MACe, MACf, r, nfT, nfM ] = GeneralandMAC( N,Nb, f, f2T,f2M, VM, VT, VMMM, VTTT, KM, KT, MM, Ma, dT, dM, VMMs, VMMns );

%% MSF
[ msfT ] = MSF( VTTs(:,s), VTTns(:,s), bo, Nb );
[ msfM ] = MSF( VMMs(:,s), VMMns(:,s), bo, Nb );
[ msfTa ] = MSF( VTTs(:,s), VTTns(:,s), 0, Nb );%all elements
[ msfMa ] = MSF( VMMs(:,s), VMMns(:,s), 0, Nb );%all elements