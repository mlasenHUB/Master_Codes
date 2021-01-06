function [ msf, f2 ] = MC1(N, e, w, s, ra, DAMP, b, g )
% Montecarlo study for the msf and max suppression for different multi
%forces of suppression
%% Creating an excitation in any/various DOF.(Just for MIMO model)
%Freqency domain:Total
for i=1:2*N
    %<alt 60 >alt 62
    if i==e %&& i<=6 %|| i==2% && i<=36 || i>40 && i<=44
       f(i)=1;
    else
       f(i)=0;
    end
end

%% K-STIFFNESS MATRIX

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

%% M-Mass Matrix
md=30*ones(N,1);%masses of the disks
mb=1*ones(N,1);%masses of the blades
di=[md;mb];
Ma=diag(di);%MASS MATRIX


%% EIGEN PROBLEM
[V,D] = eig(K,Ma);%V:eigen vectors(normalised), D:eigen values.
d=diag(D);%just the (squared) frequencies

%% Multiple RANDOM Forces

[ we,f2, Vmin ] = Suppressfindforce2( 0, f, V,s );

%% Frequency Domain
fdom=1:0.001:2*(max(D(:))^0.5);
l=length(fdom);

%% DISPLACEMENT OF MODES WITH PIEZO SUPPRESSION w/o Damping
%Find the response for all the possible positions
%WITH SUPPRESSION (f2=Piezo Force =!0)
xr=cell(1, 2*N);
xrn=cell(1, 2*N);
xtr=zeros(l, 2*N);
xtrn=zeros(l, 2*N);
phr=zeros(l, 2*N);
phrn=zeros(l, 2*N);

if DAMP==0
  
for i=1:2*N
    [ x, xt, ph ] = NoDamp( i, V,d,  f, f2, w, ra,s, fdom );
    xr{i}=x;
    xtr(:,i)=xt;  %WITH SUPPRESSION  (f2=!0) 
    phr(:,i)=ph;
end
%WITHOUT SUPPRESSION (f2=Piezo Force =0)
for i=1:2*N
    [ x, xt, ph ] = NoDamp( i, V,d,  f, 0, w, ra,s, fdom );
    xrn{i}=x;
    xtrn(:,i)=xt;  %WITHOUT SUPPRESSION (f2=0)
    phrn(:,i)=ph;
    
end

end

if DAMP==1
    
for i=1:2*N
    [ x, xt, ph ] = ModalDamp( b,g, i, V,d,  f, f2,w, ra, s, fdom );
    xr{i}=x;
    xtr(:,i)=xt;  %WITH SUPPRESSION  (f2=!0) 
    phr(:,i)=ph;
end
    
%WITHOUT SUPPRESSION
for i=1:2*N
    [ x, xt, ph ] = ModalDamp( b,g, i, V,d,  f, 0,w, ra, s, fdom );
    xrn{i}=x;
    xtrn(:,i)=xt;  %WITHOUT SUPPRESSION  (f2=0)  
    phrn(:,i)=ph;
end
    
end

%% Create vectors of the modes with and without the suppression.
for i=1:2*N

    if DAMP==0
        [val, indf]=min(abs((d(i)^0.5-fdom)));% finding the closest position in fdom to the natural frequencies
    end
    if DAMP==1
        %Imaginary:
        nr=b+g./d;
        lambda2=d.*(complex(1,nr));
        lambda=sqrt(lambda2);%Complex number
        abslambda=abs(lambda);
        [val, indf]=min(abs((abslambda(i)-fdom)));
        indfv(i)=indf;
             
    end
    

vs=xtr(indf,:); %vs contains all the displacements of the nodes at one resonance frequency WITH SUPRESSION.
vns=xtrn(indf,:); %vns contains all the displacements of the nodes at one resonance frequency WITH NO SUPRESSION.

Vs(:,i)=vs;%Vs contains all the displacementes of the nodes at all resonant frequencies WITH SUPRESSION
Vns(:,i)=vns;%Vns contains all the displacementes of the nodes at all resonant frequencies WITH NO SUPRESSION

phs=phr(indf,:); %vs contains all the displacements of the nodes at one resonance frequency WITH SUPRESSION.
phns=phrn(indf,:); %vns contains all the displacements of the nodes at one resonance frequency WITH NO SUPRESSION.

PHs(:,i)=phs;%Vs contains all the displacementes of the nodes at all resonant frequencies WITH SUPRESSION
PHns(:,i)=phns;%Vns contains all the displacementes of the nodes at all resonant frequencies WITH NO SUPRESSION
end
[ Vs, Vns, y2, y3 ] = suppeffect(Vs, Vns, DAMP, PHs, PHns );
[msf]=MSF(Vs, Vns);

%% Just the values of the suppressed mode
 msf=msf(s);
 f2=f2';
 end

