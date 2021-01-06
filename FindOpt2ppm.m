%% Find and save all the optimums for tuned and mistuned system for a model with 2DOFs per blade
%and  Piezo per mode suppression


mi=5;%4% of alternating mistuning
modes=12;
A=0; %1 if Appoach is minimum force, 0 if 1%-->MSF

%% General Control Parameters:
N=4;% number of disks
Nb=2*N;%Number of blades per disk
c=4;%Number of columns of mode shapes to display.
e=12;% Excitation Point
sine=0;%Sinusoidal excitation IN THE BLADES
sinp=0;%period of the sinusoidal exitation (ND)
DOUBLE=1; % DOUBLE forcing activated(1: yes, 0:no)
ppm=2;% piezos per mode to suppressed (1=1 piezo per mode, 2=2 piezos per mode)
w=1; %1 for windowing, else:0

ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=0;%1 for plotting, 0 for no plotting
t=0; %no torsional modes in 2nd model
p=1;%percentage of the excitation force chosen to be exerted by the piezo
loadopt=0; %(1 if yes)load the optimum from a previously saved mode
td= 0; %tip to disk piezo: 1=operating, 0=not in operation
bo=1;%1: blades only activated, 0: deactivated

gy3T=zeros(N+Nb,2*modes);%all displacement for tuned and mistuned sytesm
gy3M=zeros(N+Nb,2*modes);

forcesT=zeros(N+Nb,modes);%forces tuned
forcesM=zeros(N+Nb,modes);

for s=1:modes
    s
    %%%%%%%%%%%%%%%%
    %s= Mode to be supressed with the window(up to Nb+N)
    file=(['OPTs', num2str(s),'-ppm=2real'])';
    %%%%%%%%%%%%%%%%
    %% Mistuning Parameters
    
        
        MISTK=0; % add stiffness 0=No, 1=yes
        pk= 1*mi;%percentage of change of K (+ve add, -ve substract)
        altk=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
        dofk= nan; %degree of freedom to be changed on K (N left, right blades)
        leftbk=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
        sinmk=0;%sinusoisal mistuning
        sinpk=0;%period of the sinusoidal mistuning
        MISTM=1; %add Mass mistuning 0=No, 1=yes
        pm=1*mi; %percentage of change of Ma (+ve add, -ve substract)
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
        %%
        %TYPE A: GIVEN A FORCE , HOW MUCH I CAN SUPPRESS
        %Tuned
        nrT=b+g./dT;
        %Approach 3: 2 PIEZOS/2DOF
        %Tuned
        [ MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T, fA3T ] = SuppressgivenforceA3( VT,nrT,p,dT,f,N,Nb, td);
        nrM=b+g./dM;
        %Mistuned
        [ MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M , fA3M ] = SuppressgivenforceA3( VM,nrM,p,dM,f,N,Nb, td);
        
        %%
        %Approach 3: 2 PIEZO/2 DOF
        
        %Tuned
        [MminB3T, MnomB3T] = Suppressfinddoubleforce2( f, VT, Nb/N, td );
        %Mistuned
        [MminB3M, MnomB3M] = Suppressfinddoubleforce2( f, VM, Nb/N, td );
        %% Optimization  for MINIMUM FORCE TYPE  B
        %Approach 3
        [ MmT, MmT2, PtT2,PT ] = OptPiezo2( MminB3T );
        [ MmM, MmM2, PtM2,PM ] = OptPiezo2( MminB3M );
        idloc=nan;
        %% Type A optimization: Minimum MSF
        %Approach 3
        [ PPT, MSFfT, idlocT1A3, idlocT2A3, OPTT , DsupT] = OptPiezoA3(MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T,s, N,bo,Nb);
        [ PPM, MSFfM, idlocM1A3, idlocM2A3, OPTM , DsupM] = OptPiezoA3(MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M,s, N,bo,Nb);
        
        
        %% Creation of the Optimum suppression force 
        if A==1
            [ f2T ] = CreateSuppForce(MminB2T, MminB2M , 0, 0,idlocT, Pt, N,Nb,s,ppm);%create suppressing force.
            [ f2M ] = CreateSuppForce(MminB2T, MminB2M , MISTK, MISTM,idlocM, Pt, N,Nb,s,ppm);%create suppressing force.
        end
        
        if A==0
        f2T=PPT(:,s)-f';
        f2M=PPM(:,s)-f';
        
        f2T=f2T';
        f2M=f2M';
        end
        
        %% MODE SHAPES
        VdM=VM(1:N,:);%Disks Modal Shapes with mistuning
        VbM=VM(N+1:end,:);%Blades Modal Shapes with mistuning
        
        %Modes without mistuning
        VdT=VT(1:N,:);%Disks Modal Shapes without Mistuning(Tuned)
        VbT=VT(N+1:end,:);
        
        
        %% Frequency Domain
        fdommacx=max(sqrt(dT(N+Nb)),sqrt(dM(N+Nb)));
        fdom=0.001:0.001:1.2*(fdommacx);%Frequency domain
        
        %% Displacement of modes with and w/o mistuning
        
        
        [ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g );
        
        [ VTTT, VMMM, y2, y3 ] = suppeffect(VTT, VMM, DAMP, PHT, PHM );%here I used suppeffect function(which was created for picking up the displacement from the FRF for suppression comparison)
        %% TYPE B optimization: Suppression effect on the mistuned model
        
        %tuned
        [ VTs, VTns, PHTs, PHTns, abslambdaT, xrTs, xtrTs, phrTs, xrTns, xtrTns, phrTns ] = SuppComparison( 1,dT, fdom, VT, f, f2T, w, ra, s, b, g );%Damping set to 1(activated)
        
        [ VTTs, VTTns, y2sT, y3sT ] = suppeffect(VTs, VTns, 1, PHTs, PHTns );
        
        %mistuned
        [ VMs, VMns, PHMs, PHMns, abslambdaM, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2M, w, ra, s, b, g );%Damping set to 1(activated)
        
        [ VMMs, VMMns, y2sM, y3sM ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );
        
        %% TYPE A Optimisation: Suppression effect on the mistuned model
        [ r ] = VMMs./VMMns;
        [ DTT, msfFRFT, msfCFT,oedT, oebmT,oebT, rCFT , fsT ] = suppeffectA3(s, idlocT1A3, idlocT2A3,f, f2T, MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, VTTs,PPT, OPTT,r, 'TUNED',bo, Nb );
        [ DMM, msfFRFM, msfCFM,oedM, oebmM,oebM, rCFM , fsM ] = suppeffectA3(s, idlocM1A3, idlocM2A3,f, f2M, MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, VMMs,PPM, OPTM,r, 'MISTUNED',bo, Nb );%Mistuned
       
        
        %% Usefull vectors    
        gy3T(:, (2*(s-1)+1):(2*s))=y3sT(:,(2*(s-1)+1):(2*s));
        gy3M(:, (2*(s-1)+1):(2*s))=y3sM(:,(2*(s-1)+1):(2*s));
        
        forcesT(:,s)=fsT;%forces tuned
        forcesM(:,s)=fsM;%forces tuned
        
        
        
        save(file)
end

%%
        msfT=OPTT(:,3); %all optimium msf per mode
        msfM=OPTM(:,3);