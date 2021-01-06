%% Model 2: 1 Mass disk + 2 masses blades in series: can compare the TUNED cases for Linear and Torsional Modeling.
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
t=0; %0=no torsion, 1= torsion with 2 masses in parallel, 2=torsion with 2 masses in series.
%% K-STIFFNESS MATRIX (LINEAR AND TORSION)

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kbl=1000*ones(N,1);%stiffness of the blade 'i' LEFT
kbr=1000*ones(N,1);%stiffness of the blade 'i' RIGHT
kg=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'

kdrt=1000000*ones(N,1);%stiffness to the right spring of disk 'i' TORSION
kdlt=fliplr(kdr);%stiffness to the left spring of disk 'i' TORSION
kblt=100*ones(N,1);%stiffness of the blade 'i' LEFT in TORSION
kbrt=100*ones(N,1);%stiffness of the blade 'i' RIGHT in TORSION
kgt=10000*ones(N,1);%stiffness of the disk-shaft spring of disk 'i'


H= diag(kdl+kdr+kbl+kg);%elements of K from the equations of motion of the disks
H2= diag(-kdr(1:N-1),1);%concatenate elements from equations of motion of the blades
H=H+H2;
H(1,N)=-kdr(N);
H=H+triu(H,1)';

%torsion
Ht= diag(kdlt+kdrt+kblt+kgt);%elements of K from the equations of motion of the disks
H2t= diag(-kdrt(1:N-1),1);%concatenate elements from equations of motion of the blades
Ht=Ht+H2t;
Ht(1,N)=-kdrt(N);
Ht=Ht+triu(Ht,1)';

D0=diag(zeros(1,N));
D1=diag(-kbl);
D2=diag(kbl+kbr);
D3=diag(-kbr);
D4=diag(kbr);

D1t=diag(-kblt);
D2t=diag(kblt+kbrt);
D3t=diag(-kbrt);
D4t=diag(kbrt);

K= [H D1 D0; D1 D2 D3; D0 D3 D4]; %linear
Kt= [Ht D1t D0; D1t D2t D3t; D0 D3t D4t]; %Torsion

%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mbl=1*ones(N,1);%masses of the blades LEFT
mbr=1*ones(N,1);%masses of the blades RIGHT

rad=2;%radious of the plate representing the blades.
id=md*((1000*rad)^2/4);%masses of the disks
ibl=mbl*(rad^2/4);%inertia of the LEFT plate
ibr=mbr*(rad^2/4);%inertia of the RIGHT plate


Ma=diag([md;mbl;mbr]);%MASS MATRIX
Ia=diag([id; ibl; ibr]);%INERTIA MATRIX
%MI=diag([md;mbl;mbr;ibl;ibr]); %MASS & INERTIA MATRIX
%% EIGEN PROBLEM
[VT,DT] = eig(K,Ma);%With Tuned linear
dT=diag(DT);
[VTt,DTt] = eig(Kt,Ia);%With Tuned torion just blades
dTt=diag(DTt);

%% MODE SHAPES
%Tuned Modes
VdT=VT(1:N,:);%Disks 
VbT=VT(N+1:end,:);%LEFT and RIGHT blades LINEAR

VdTt=VTt(1:N,:);%Disks 
VbTt=VTt(N+1:end,:);%LEFT and RIGHT blades TORSION


%% PLOTTING Mode Shapes and Natural Frequencies vs Nodal Diameters
if plotmodes==1
Plotmodes(c, dT, VT,VdT, VbT, t )%Linear
Plotmodes(c, dTt, VTt,VdTt, VbTt, t )%Torsion
end
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
[Mmin2, Mnom2] = Suppressfindforce( f, VT );

Mmin2c=Mmin2.*(abs(Mmin2)<1.001);% Mmin2 constrain to forces less than the excitation(unity, 0.001 error)
P= Mmin2(pp,s);%See matrix Mmin2% Force applied
for i=1:Nb+N
    
    if  i==pp% && i<=36 || i>40 && i<=44
        f2(i)=P;
    else
        f2(i)=0;
    end
end

%% Responses to Forcing with and w/o damping
fdom=0.001:0.001:2*(max(DT(:))^0.5);%Frequency domain


[ x, xtu, phu ] = NoDamp( r, VT,dT,  f, f2, w, ra,s, fdom );%undamped
[ xd, xtd, phd ] = ModalDamp( b,g, r, VT,dT,  f, f2,w, ra, s, fdom );%damped

PlotTF( xd,xtd, f2, pp, fdom, r )%Plot transference function and its modes contributions.

PlotMain(x,xtu,xtd, phu,phd,dT, fdom, VdT, VbT, f, f2)%Plot the summary of results: Total transference function, peaks
%% DISPLACEMENT OF MODES WITH PIEZO SUPPRESSION with and w/o Damping
[ Vs, Vns, PHs, PHns, abslambda, xr, xtr, phr, xrn, xtrn, phrn ] = SuppComparison( DAMP, dT , fdom, VT, f, f2, w, ra, s, b, g );

PlotTFpoints( c, xr, xtr,phr, fdom, 'WITH SUPPRESSION')  %Suppressed 
PlotTFpoints( c, xrn, xtrn, phrn, fdom, 'WITHOUT SUPPRESSION')  % Not suppressed

[ Vs, Vns, y2, y3 ] = suppeffect(Vs, Vns, DAMP, PHs, PHns );
barcomparison( c, Vs, Vns, 'Suppressed', 'Unsuppressed');
%% Suppression
[msf]=MSF(Vs, Vns);
figure
[ Rd ] = ratio( Vs, Vns );
