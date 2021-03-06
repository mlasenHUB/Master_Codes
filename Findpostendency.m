%% FIND Tendency: finds if the piezos tend to be to higher positions as the frequency increase

%trying 10 cases of mistuning from 1 to 10%
mistcases=1;
OPT=zeros(12,4,mistcases); %opt positions and msf for every level of mistuning
OPM=zeros(12,4,mistcases);
modes=1; %modes to study


    
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

ra=0.4; %Range of Piezo Forcing Frequency
DAMP=1; %Damping model, 0:No Damping, 1: Modal Damping(C=H=bK+gM)
b=0.0005;
g=0.08;
MULT=0;%Multiple forcing yes=1, no=0;
plotmodes=0;%1 for plotting, 0 for no plotting
t=0; %no torsional modes in 2nd model
p=1;%percentage of the excitation force chosen to be exerted by the piezo
loadopt=0; %(1 if yes)load the optimum from a previously saved mode
td= 0; %tip to disk piezo: 1=operating, 0=not in operation
bo=1; %analise blades-only in the Modal Scale Factor (1:yes, 0:no)


for s=1:modes
    %%%%%%%%%%%%%%%%
    %s= Mode to be supressed with the window(up to Nb+N)
    file=(['Ts', num2str(s)])';%2:blades only(m-2-s, k-2-s) 1:all dofs
    %%%%%%%%%%%%%%%%
    %% Mistuning Parameters
    for mi=1:mistcases
        s,mi
        MISTK=0; % add stiffness 0=No, 1=yes
        pk= 1*mi;%percentage of change of K (+ve add, -ve substract)
        altk=1;%1: activate, 0, deactivated(alternate mistuning +1, -1, etc.)
        dofk= nan; %degree of freedom to be changed on K (N left, right blades)
        leftbk=0;%1 if the mistuning is applied to the LEFT mass blade (0 for the RIGHT)
        sinmk=0;%sinusoisal mistuning
        sinpk=0;%period of the sinusoidal mistuning
        MISTM=0; %add Mass mistuning 0=No, 1=yes
        pm=0*mi; %percentage of change of Ma (+ve add, -ve substract)
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
        f(e)=1;f(9)=1;f(10)=1;f(11)=1;
        %%
        %TYPE A: GIVEN A FORCE , HOW MUCH I CAN SUPPRESS
        %Approach 3: 2 PIEZOS/2DOF
        %Tuned
        nrT=b+g./dT;
        [ MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T, fA3 ] = SuppressgivenforceA3( VT,nrT,p,dT,f,N,Nb, td);
        %Mistuned
        nrM=b+g./dM;
        [ MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M , fA3 ] = SuppressgivenforceA3( VM,nrM,p,dM,f,N,Nb, td);
        
        
        
        %% Type A optimization: Minimum MSF
        %Approach 3
        
        [ PPT, MSFfT, idlocT1A3, idlocT2A3, OPTT ] = OptPiezoA3(MA3T, MppA3T, MnpA3T, MnnA3T,MpnA3T, FppA3T, FnpA3T,FnnA3T, FpnA3T,s, N,bo,Nb);
        [ PPM, MSFfM, idlocM1A3, idlocM2A3, OPTM ] = OptPiezoA3(MA3M, MppA3M, MnpA3M, MnnA3M,MpnA3M, FppA3M, FnpA3M,FnnA3M, FpnA3M,s, N,bo,Nb);
        
        OPT(:,:,mi)=OPTT;
        OPM(:,:,mi)=OPTM;
        
        
    end
    
    save(file)
end