%% Model 2: 1 Mass disk + 2 masses blades.
clear all
close all
%% General Control Parameters:
N=4;% number of disks
Nb=2*N;%Number of blades per disk (in total N or 2*N)
c=4;%Number of columns of mode shapes to display.
e=12;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=2;%period of the sinusoidal exitation (ND)
%r=7;% Response Point/Drive Point
DOUBLE=1; % DOUBLE forcing activated(1: yes, 0:no)
ppm=2;% piezos per mode to suppressed (1=1 piezo per mode, 2=2 piezos per mode)
w=1; %1 for windowing, else:0
s=8; %Mode to be supressed with the window(up to Nb+N)
ra=0.01; %Half Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;%0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=0;%1 for plotting, 0 for no plotting
plotdeandmc=0;%1 for plotting damping effect and modal contribution(0 for no)
t=0; %no torsional modes in 2nd model
p=1;%percentage of the excitation force chosen to be exerted by the piezo
loadopt=0; %(1 if yes)load the optimum from a previously saved mode
td= 0; %tip to disk piezo: 1=operating, 0=not in operation
bo=1;%analise blades-only in the Modal Scale Factor (1:yes, 0:no)
msfforce=1; %check the effect in the suppression with the 1% of excitation force in the FRF (1 yes, 0 no)
%% Mistuning Parameters

MISTK=0; % add stiffness 0=No, 1=yes 
pk= 1;%percentage of change of K (+ve add, -ve substract)
altk=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
dofk= nan; %degree of freedom to be changed on K (N left, right blades)
leftbk=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinmk=0;%sinusoisal mistuning
sinpk=0;%period of the sinusoidal mistuning

MISTM=0; %add Mass mistuning 0=No, 1=yes
pm=5; %percentage of change of Ma (+ve add, -ve substract)
altm=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
dofm=0; %degree of freedom to be changed on M (Nb+N Masses)
leftbm=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
sinm=0;%sinusoisal mistuning(1 if yes, 0 if no)
sinpm=0;%period of the sinusoidal mistuning

%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kbl=1000*ones(N,1);%stiffness of the blade 'i' LEFT
kbr=1000*ones(N,1);%stiffness of the blade 'i' RIGHT
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'

if MISTK==1
    if leftbk==1
    [ kblM, kbl ] = MistuningStiffness( kbl, pk, dofk, N, sinmk, sinpk, altk  );
    kbrM=kbr;
    else
    [ kbrM, kbr ] = MistuningStiffness( kbr, pk, dofk, N, sinmk, sinpk, altk  );
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
mbl=1*ones(N,1);%masses of the blades LEFT
mbr=1*ones(N,1);%masses of the blades RIGHT
di=[md;mbl;mbr];
Ma=diag(di);%MASS MATRIX

if MISTM==1
    [ MM, Maor ] = MistuningMass( Ma, pm , dofm, N, sinm, sinpm,leftbm, altm );
else 
    MM=Ma;
end 
%% EIGEN PROBLEM
[VT,DT] = eig(KT,Ma);% Tuned
dT=diag(DT);
[VM,DM] = eig(KM,MM);%MisTuned
dM=diag(DM);
%% Creating an excitation in any/various DOF.(Just for MIMO model)

f= zeros(1, Nb+N);
f(e)=1;
f(9)=1;f(11)=1;f(10)=1;

%% Finding the best positions, forcing values to suppress certain modes

%TYPE A: GIVEN A FORCE , HOW MUCH I CAN SUPPRESS
% Approach 1 : 1 PIEZO/1 DOF

%Tuned
nrT=b+g./dT;
[ MA1T, MpA1T, MnA1T ] = Suppressgivenforce( VT,nrT,p,dT,f);
%Mistuned
nrM=b+g./dM;
[ MA1M, MpA1M, MnA1M ] = Suppressgivenforce( VM,nrM,p,dM,f);

% Approach 2: 1 PIEZO/2 DOF
%Tuned
[ MA2T, MpA2T, MnA2T, FpA2T, FnA2T,fA2T ] = SuppressgivenforceA2( VT,nrT,p,dT,f,N,Nb,td);
%Mistuned
[ MA2M, MpA2M, MnA2M,  FpA2M, FnA2M,fA2M ] = SuppressgivenforceA2( VM,nrM,p,dM,f,N,Nb,td);

%Approach 3: 2 PIEZOS/2DOF
%Tuned
[ MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T, fA3T ] = SuppressgivenforceA3( VT,nrT,p,dT,f,N,Nb, td);
%Mistuned
[ MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M , fA3M ] = SuppressgivenforceA3( VM,nrM,p,dM,f,N,Nb, td);

%TYPE B: FIND THE MINIMUM FORCE TO SUPPRESS A MODE COMPLETELY

% Approach 1 : 1 PIEZO/1 DOF
%Tuned
[MminB1T, MnomB1T] = Suppressfindforce( f, VT );
MminB1Tc=MminB1T.*(abs(MminB1T)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)

%Mistuned
[MminB1M, MnomB1M] = Suppressfindforce( f, VM );
MminB1Mc=MminB1M.*(abs(MminB1M)<1.001);% Mmin2 Constrain to forces less than the excitation(unity, 0.001 error)

%Approach 2: 1 PIEZO/2 DOF
%Tuned
[MminB2T, MnomB2T] = Suppressfinddoubleforce( f, VT, 2,1 );%td=0
%Mistuned
[MminB2M, MnomB2M] = Suppressfinddoubleforce( f, VM, 2,1 );

%Approach 3: 2 PIEZO/2 DOF

%Tuned
[MminB3T, MnomB3T] = Suppressfinddoubleforce2( f, VT, Nb/N, td );
%Mistuned
[MminB3M, MnomB3M] = Suppressfinddoubleforce2( f, VM, Nb/N, td );

%% Optimization  for MINIMUM FORCE
%T=Tuned; M= Mistuned

% TYPE B optimization: MINIMUM FORCE

%Approach 2 
if DOUBLE==1 && ppm==1
    [ FtT, PtT , MfT, idlocT ] = OptPiezo( MminB2T,s );
    [ FtM, PtM, MfM, idlocM ] = OptPiezo( MminB2M,s );

    Pt=nan;

    %PlotOpt(c,FtT, PtT, MfT,'TUNED' )
    %PlotOpt(c,FtM, PtM, MfM,'MISTUNED' )
end

%Approach 3
if DOUBLE==1 && ppm==2%TWO piezos per mode
    [ MmT, MmT2, PtT2,PT ] = OptPiezo2( MminB3T );
    [ MmM, MmM2, PtM2,PM ] = OptPiezo2( MminB3M );   
    idloc=nan;
end

%% Optimization of the LOCATION of the suppression force for minimum MSF(FORCE is CHOSEN)
%Type A optimization: Minimum MSF

%Approach 1

if DOUBLE==0
   [ PPT, MSFfT, idlocTA1 ] = OptPiezo0( f, MA1T, MpA1T, MnA1T,s, bo, Nb );%tuned
   [ PPM, MSFfM, idlocMA1] = OptPiezo0( f, MA1M, MpA1M, MnA1M,s, bo, Nb );%mistuned
end

%Approach 2

if DOUBLE==1 && ppm==1
    
    [ PPT, MSFfT, idlocTA2 ] = OptPiezoA2( MA2T, MpA2T, MnA2T, FpA2T, FnA2T, s, N, bo,Nb );
    [ PPM, MSFfM, idlocMA2 ] = OptPiezoA2( MA2M, MpA2M, MnA2M,FpA2M, FnA2M, s, N, bo,Nb );
end

%Approach 3

if DOUBLE==1 && ppm==2
    [ PPT, MSFfT, idlocT1A3, idlocT2A3, OPTT , DsupT] = OptPiezoA3(MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T,s, N,bo,Nb);
    [ PPM, MSFfM, idlocM1A3, idlocM2A3, OPTM , DsupM] = OptPiezoA3(MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M,s, N,bo,Nb);
    
end

%% Creation of the Optimum suppression force TYPE B(from the minimum)

%Just comparing the suppressing effect on the model with mistuning (neglect
%the tuned model)
if DOUBLE==0
    %Tuned
    f2T=zeros(1,N+Nb);
    [mvT, idlocT]=min(abs(MminB1T(:,w)));
    f2T(idlocT)=MminB1T(idlocT,s);%create suppressing force.
    %Mistuned
    f2M=zeros(1,N+Nb);
    [mvM, idlocM]=min(abs(MminB1M(:,w)));
    f2M(idlocM)=MminB1M(idlocM,s);%create suppressing force
end
if DOUBLE==1 && ppm==1%ONE piezo per mode
    [ f2T ] = CreateSuppForce(MminB2T, MminB2M , 0, 0,idlocT, Pt, N,Nb,s,ppm);%create suppressing force.
    [ f2M ] = CreateSuppForce(MminB2T, MminB2M , MISTK, MISTM,idlocM, Pt, N,Nb,s,ppm);%create suppressing force.
    if msfforce==1
        f2T=PPT(:,s)-f';
        f2M=PPM(:,s)-f';
        f2T=f2T';
        f2M=f2M';
    end
end
if DOUBLE==1 && ppm==2%TWO piezo per mode
    [ f2T ] = CreateSuppForce(MminB3T, MminB3M , 0, 0,idloc, PtT2, N,Nb,s,ppm);%create suppressing force.
    [ f2M ] = CreateSuppForce(MminB3T, MminB3M , MISTK, MISTM,idloc, PtM2, N,Nb,s,ppm);%create suppressing force.
    if msfforce==1
        f2T=PPT(:,s)-f';
        f2M=PPM(:,s)-f';
        f2T=f2T';
        f2M=f2M';
    end

end
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
fdommacx=max(sqrt(dT(N+Nb)),sqrt(dM(N+Nb)));
fdom=0.001:0.001:1.2*(fdommacx);%Frequency domain
%% PLOTTING DAMPING EFFECT AND MODES CONTRIBUTION AT 1 LOCATION
if plotdeandmc==0
    r=12;
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
%% TYPE B optimization: Suppression effect on the mistuned model

%tuned

[ VTs, VTns, PHTs, PHTns, abslambdaT, xrTs, xtrTs, phrTs, xrTns, xtrTns, phrTns ] = SuppComparison( 1,dT, fdom, VT, f, f2T, w, ra, s, b, g );%Damping set to 1(activated)
PlotTFpoints( c, xrTs, xtrTs, phrTs, fdom, 'TUNED-SUPPRESSED')
PlotTFpoints( c, xrTns, xtrTns, phrTns, fdom, 'TUNED-UNSUPPRESSED')

[ VTTs, VTTns, y2sT, y3sT ] = suppeffect(VTs, VTns, 1, PHTs, PHTns );
barcomparison( c, VTTs, VTTns, 'Tuned-Suppressed', 'Tuned-Unsuppressed');

%mistuned

[ VMs, VMns, PHMs, PHMns, abslambdaM, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2M, w, ra, s, b, g );%Damping set to 1(activated)
PlotTFpoints( c, xrMs, xtrMs, phrMs, fdom, 'MISTUNED-SUPPRESSED')
PlotTFpoints( c, xrMns, xtrMns, phrMns, fdom, 'MISTUNED-UNSUPPRESSED')

[ VMMs, VMMns, y2sM, y3sM ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );
barcomparison( c, VMMs, VMMns, 'Mistuned-Suppressed', 'Mistuned-Unsuppressed');

%% General results: ratio of diplacement, MACs, natural frequencies and suppression

[MACe, MACf, r, nfT, nfM ] = GeneralandMAC( N,Nb, f, f2T,f2M, VM, VT, VMMM, VTTT, KM, KT, MM, Ma, dT, dM, VMMs, VMMns );

%% TYPE A Optimisation: Suppression effect on the mistuned model 
if DOUBLE==0
    [ DTT, msfFRFT, msfCFT,oedT, oebmT,oebT, rCFT, fsT ] = suppeffect2( s, idlocTA1,f, f2T, MA1T,MpA1T,MnA1T, VTTs,PPT,r, 'Tuned',bo, Nb );
    [ DMM, msfFRFM, msfCFM,oedM, oebmM,oebM, rCFM, fsM ] = suppeffect2( s, idlocMA1,f, f2M, MA1M,MpA1M,MnA1M, VMMs,PPM,r, 'Mistuned',bo, Nb );
end
if DOUBLE==1 && ppm==1
    [ DTT, msfFRFT, msfCFT,oedT, oebmT,oebT, rCFT, fsT ] = suppeffect2( s, idlocTA2,f, f2T, MA2T,MpA2T,MnA2T, VTTs,PPT,r, 'Tuned',bo, Nb );
    [ DMM, msfFRFM, msfCFM,oedM, oebmM,oebM, rCFM, fsM ] = suppeffect2( s, idlocMA2,f, f2M, MA2M,MpA2M,MnA2M, VMMs,PPM,r, 'Mistuned',bo, Nb );
end
if DOUBLE==1 && ppm==2
    [ DTT, msfFRFT, msfCFT,oedT, oebmT,oebT, rCFT , fsT ] = suppeffectA3(s, idlocT1A3, idlocT2A3,f, f2T, MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, VTTs,PPT, OPTT,r, 'TUNED',bo, Nb );
    [ DMM, msfFRFM, msfCFM,oedM, oebmM,oebM, rCFM , fsM ] = suppeffectA3(s, idlocM1A3, idlocM2A3,f, f2M, MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, VMMs,PPM, OPTM,r, 'MISTUNED',bo, Nb );%Mistuned
end

%% Retrieve the optimium forces for another case to search its effect on the current mode
if loadopt==1
    [ DpT, sor , msfpT, msfcT, bodpT, bodcT ] = fixdpiezoloc( 'OPTs6boT', 0, VT, nrT, dT,s, f,DTT, OPTT, PPT, DOUBLE, ppm, 'TUNED' , bo, Nb);
    [ DpM, sor , msfpM, msfcM, bodpM, bodcM ] = fixdpiezoloc( 'OPTs6boM', 1, VM, nrM, dM,s, f,DMM, OPTM, PPM, DOUBLE, ppm, 'MISTUNED', bo, Nb );
end

%% Usefull vectors
msfT=max(MSFfT);
msfM=max(MSFfM);

for i=1:12
    msfTt(i)=max(msfT(1,:,i));
    msfMt(i)=max(msfM(1,:,i));
end
figure
plot(msfTt)
hold on
plot(msfMt)
grid on
