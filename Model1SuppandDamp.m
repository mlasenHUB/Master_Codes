clear all
close all
%% General Control Parameters:
N=4;%2*N mode shapes.
c=4;%Number of columns of mode shapes to display.
e=6;% Excitation Point
r=7;% Response Point/Drive Point
pp=6;%Piezo position
w=1; %1 for windowing, else:0
s=8; %Mode to be supressed with the window(up to 2*N
ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=1;%1 for plotting, 0 for no plotting
t=0;%no torsion modes
%% Mistuning Parameters

MISTK=0; % add stiffness 0=No, 1=yes 
pk= 70;%percentage of change of K (+ve add, -ve substract)
dofk= 1; %degree of freedom to be changed on K (N blades)

MISTM=0; %add Mass mistuning 0=No, 1=yes
pm=70; %percentage of change of Ma (+ve add, -ve substract)
dofm=9; %degree of freedom to be changed on M (2N Masses, start from N+1)

%% Creating an excitation in any/various DOF.(Just for MIMO model)
%Freqency domain:Total
for i=1:2*N
    %<alt 60 >alt 62
    if i==e %|| i==e+2% && i<=36 || i>40 && i<=44
       f(i)=1;
    
    else
       f(i)=0;
    end
end

 %      f(6)=-1;f(8)=-1;

%% K-STIFFNESS MATRIX

kdr=10000*ones(N,1);%stiffness to the right spring of disk 'i'
kdl=fliplr(kdr);%stiffness to the left spring of disk 'i'
kb=1000*ones(N,1);%stiffness of the blade 'i'
%Add mistuning in a blade
if MISTK==1
    [ kb, kbor ] = MistuningStiffness( kb, pk, dofk, N, 0, 1  );
else
    kbor=kb;
end

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

%K without mistuning

d= kdl+kdr+kbor+kg;%elements of K from the equations of motion of the disks
d= [d;kbor];%concatenate elements from equations of motion of the blades

KD= diag(d);%create basic Stiffness matrix diagonal, without coupling elements 
Kr=diag([(-kdr(1:N-1));z],1);%to the right of the diagonal
Krc=diag((-kbor), N);%to the right corner of the diagonal
Kor=Kr+Krc;
Kor(1,N)=-kdr(N);%term due to the left of the first disk
Kor=Kor+Kor.';%Symmetric
Kor=Kor+KD;%STIFFNESS MATRIX 


%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mb=1*ones(N,1);%masses of the blades
di=[md;mb];
Ma=diag(di);%MASS MATRIX

if MISTM==1
    [ Ma, Maor ] = MistuningMass( Ma, pm , dofm, N, 0, 1 );
else 
    Maor=Ma;
end  

%% EIGEN PROBLEM
[V,D] = eig(K,Ma);%V:eigen vectors(normalised), D:eigen values.

[Vor, Dor]=eig(Kor, Maor); %eigen problem with the mistuned case
d=diag(D);%just the (squared) frequencies
dor=diag(Dor);

%{
%THRESHOLD AT NODAL DIAMETER
for i=1:2*N
    
    maxV=max(abs(V(:,i)));
    Vw(:,i)=V(:,i).*(abs(V(:,i))> maxV/2);
end
[Vw,D] = eig(K,Ma);%V With filtering non-zero values at nodal diameters.
%}

%% Check Normalisation
Identity=V'*(Ma*V);
ResFreqs=V'*(K*V);
%% Finding the best positions, forcing values to suppress certain modes

%First Approach : Mmin: find how much I can suppress a mode with a given force of magnitude
%F.
F=1;
[Mmin]=Suppressgivenforce( c,f, Vor,F );

%Second Approach : Mmin2: What is the best force magnitude, to place at any position, to suppress a certain mode
[Mmin2, Mnom2] = Suppressfindforce( f, Vor );

Mmin2c=Mmin2.*(abs(Mmin2)<1.001);% Mmin2 constrain to forces less than the excitation(unity, 0.001 error)


P= Mmin2(pp,s);%See matrix Mmin2% Force applied
for i=1:2*N
    
    if  i==pp% && i<=36 || i>40 && i<=44
        f2(i)=P;
    else
        f2(i)=0;
    end
end
%Multiple forcing for suppression

[ we,fs, Vmin ] = Suppressfindforce2( 0, f, V,s );
if MULT==1
    f2=fs;
end
%% MODE SHAPES
Vnd=V(1:N,:);%Disks Modal Shapes with mistuning
Vnb=V(N+1:2*N,:);%Blades Modal Shapes with mistuning

%Modes without mistuning
VndWOM=Vor(1:N,:);%Disks Modal Shapes WITH OUT Mistuning
VnbWOM=Vor(N+1:2*N,:);


%% PLOTTING Mode Shapes and Natural Frequencies vs Nodal Diameters
if plotmodes==1
Plotmodes(c, d, V,Vnd, Vnb, t )
Plotmodes(c, dor, Vor,VndWOM, VnbWOM, t )
end
figure
[ nfd ] = natfreqND( dor, N, N);
%% Responses to Forcing with and w/o damping
fdom=0.001:0.001:2*(max(D(:))^0.5);%Frequency domain


[ x, xtu, phu ] = NoDamp( r, V,d,  f, f2, w, ra,s, fdom );%undamped
[ xd, xtd, phd ] = ModalDamp( b,g, r, V,d,  f, f2,w, ra, s, fdom );%damped

PlotTF( xd,xtd, f2, pp, fdom, r )%Plot transference function and its modes contributions.

PlotMain(x,xtu,xtd, phu,phd,d, fdom, Vnd, Vnb, f, f2)%Plot the summary of results: Total transference function, peaks modes, forces(excitations and piezo), and damping comparison
%% DISPLACEMENT OF MODES WITH PIEZO SUPPRESSION with and w/o Damping
[ Vs, Vns, PHs, PHns, abslambda, xr, xtr, phr, xrn, xtrn, phrn ] = SuppComparison( DAMP, d , fdom, V, f, f2, w, ra, s, b, g );


PlotTFpoints( c, xr, xtr,phr, fdom, 'WITH SUPPRESSION')  %Suppressed 
PlotTFpoints( c, xrn, xtrn, phrn, fdom, 'WITHOUT SUPPRESSION')  % Not suppressed

[ Vs, Vns, y2, y3 ] = suppeffect(Vs, Vns, DAMP, PHs, PHns );
barcomparison( c, Vs, Vns, 'Suppressed', 'Unsuppressed');

%[ SF , SFv] = Suppression( R, s, N );

%% 
%Plotting Mode contribution with and without piezo forcing
Vsd=Vs(1:N,:);%Disks Modal Shapes
Vsb=Vs(N+1:2*N,:);%Blades Modal Shapes

Vnsd=Vns(1:N,:);%Disks Modal Shapes
Vnsb=Vns(N+1:2*N,:);%Blades Modal Shapes

%PLOTTING Mode Shapes
%Plotmodes( d, Vs,Vsd, Vsb )%Suppressed modes
%Plotmodes( d, Vns,Vnsd, Vnsb )%NO Suppressed modes

%% SUPPRESSION FACTOR
%Local Supression Factor: SFl
%SFl=(100/length(Vsp))*norm(Vsp-Vs)/norm(Vs);%w.r.t exited
%SFl2= (100/length(Vsp))*norm((Vsp-V(:,s)))/norm(V(:,s));%w.r.t mode shape only excited (filtered)
%SFl3= (100/length(Vsp))*norm((Vsp-Vw(:,s)))/norm(Vw(:,s));%w.r.t mode shape only excited (filtered)
[msf]=MSF(Vs, Vns);
figure
[ Rd ] = ratio( Vs, Vns );
%%
%Displaying values of interest
fprintf('Disk Mass : %g \n', md(1));
fprintf('Blade Mass : %g \n', mb(1));
fprintf('Response Point : %g \n',r);
fprintf('Excitation Point : %g \n',e);
fprintf('Piezo position : %g \n',pp);

 
