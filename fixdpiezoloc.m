function [ Dp, sor, minmsf, msfc, bodp, bodc ] = fixdpiezoloc( file, t, V, nr, d,s, f,D,OPTc, PP, DOUBLE, ppm , label, bo, Nb)
%fixdpiezoloc finds the optimal force to suppress a mode with the position
%of the piezos fixed for another mode optimisation.
%   file= name of the file to recall the optimal position(s) from.
%   t= tuning: 0=tuned, 1=mistuned
%   V=eigenvectors
%   nr=damping factor
%   d= eigenvalues
%   f=excitation force

if t==0
    load(file, 'f2T')
    
    OPTstruct=load(file, 'OPTT');
    OPTcell=struct2cell(OPTstruct);
    OPT=cell2mat(OPTcell);
    
end

if t==1
    load(file, 'f2M')
    
    OPTstruct=load(file, 'OPTM');
    OPTcell=struct2cell(OPTstruct);
    OPT=cell2mat(OPTcell);
end

sorstruct=load(file, 's');
sorcell=struct2cell(sorstruct);
sor=cell2mat(sorcell);

fAstruct=load(file, 'fA3');
fAcell=struct2cell(fAstruct);
fA=cell2mat(fAcell);


%OPTIMUM LOCATIONS FROM PREVIOUS ITERATION
idloc1=OPT(sor,1);%first piezo
idloc2=OPT(sor,2);%second piezo

n=length(f);
%generate all possible forces executable in those positions
if DOUBLE ==1 && ppm==2
    fp1=fA(idloc1,:); %force 1st piezo
    fp2=fA(idloc2,:); %force 2nd piezo
    %all possible combination of these 2 forces
    fppp=f+fp1+fp2;
    fpnp=f+fp1-fp2;
    fpnn=f-(fp1-fp2);
    fppn=f-(fp1+fp2);
    
      
    %displacement caused by each force at each DOF
    a=f*V;
    app=fppp*V;
    anp=fpnp*V;
    ann=fpnn*V;
    apn=fppn*V;
    
    
    den=1i*nr(s)*d(s);
    M=zeros(n,1);
    Mpp=zeros(n,1);
    Mnp=zeros(n,1);
    Mnn=zeros(n,1);
    Mpn=zeros(n,1);
    for ii=1:n
        for m=1:n
            M(ii)=M(ii)+ (V(ii,m)*a(m))/(d(m)-d(s)+den);   %wothout suppression force
            Mpp(ii)= Mpp(ii)+ (V(ii,m)*app(m))/(d(m)-d(s)+den);%k: mode, ii:force position.(a(k):'pseudo' modal constant
            Mnp(ii)= Mnp(ii)+ (V(ii,m)*anp(m))/(d(m)-d(s)+den);
            Mnn(ii)= Mnn(ii)+ (V(ii,m)*ann(m))/(d(m)-d(s)+den);
            Mpn(ii)= Mpn(ii)+ (V(ii,m)*apn(m))/(d(m)-d(s)+den);
        end
    end
    
    Ma=abs(M).*sign(sin(deg2rad(angle(M))));
    Mpp=abs(Mpp).*sign(sin(deg2rad(angle(Mpp))));
    Mnp=abs(Mnp).*sign(sin(deg2rad(angle(Mnp))));
    Mnn=abs(Mnn).*sign(sin(deg2rad(angle(Mnn))));
    Mpn=abs(Mpn).*sign(sin(deg2rad(angle(Mpn))));
    
    %Optimization: minimum MSF
    
    msfpp=MSF(Mpp, Ma, bo, Nb);%1
    msfnp=MSF(Mnp, Ma, bo, Nb);%2
    msfnn=MSF(Mnn, Ma, bo, Nb);%3
    msfpn=MSF(Mpn, Ma, bo, Nb);%4
    
    msf=[msfpp, msfnp, msfnn, msfpn];
    
    [minmsf,pos]=min(msf);
        
    if pos==1
        f2s=fppp-f;
        DsCFp= Mpp; %displacement of the suppressed modes with th Chosen force from Previous iteration
    end
    if pos==2
        f2s=fpnp-f;
        DsCFp= Mnp; %displacement of the suppressed modes with th Chosen force from Previous iteration
    end
    if pos==3
        f2s=fpnn-f;
        DsCFp= Mnn; %displacement of the suppressed modes with th Chosen force from Previous iteration
    end
    if pos==4
        f2s=fppn-f;
        DsCFp= Mpn; %displacement of the suppressed modes with th Chosen force from Previous iteration
    end
    Dus= D(:,1); %unsuppressed displacement
    DsCF=D(:,3); %suppresses with optimum force for current optimum
    
    Dp=[Dus , DsCF , DsCFp];%Displacements comparison of the optimum from the current mode with the optimum from a previous iteration
end

msfc=OPTc(s,3); %minimum msf in current iteration.
bodc=DsCF./Dus;
bodp=DsCFp./Dus;

%{
figure('units','normalized','outerposition',[0 0 1 1])
suptitle(label)
subplot(3,1,1)
    bar(D)
    xlabel('DOF')
    ylabel('Displacement')
    legend('UNSUPPRESED Displacement','SUPPRESSED WITH OPTIMUM CHOSE FORCE', 'SUPPRESSED WITH OPTIMUM CHOSEN FORCE FOR GIVEN PIEZO POSITIONS')
    title([' Current Mode: ', num2str(s),'- Original Mode: ', num2str(sor),'-- MSF -current: ',num2str(msfc),'; MSF - original ', num2str(minmsf)])

fs=PP(:,s)-f';
subplot(3,1,2)
    stem(fs)
    legend('Optimum force with chosen magnitude- current optimisation')
    xlabel('Element')
   
subplot(3,1,3)

    stem(f2s)    
    legend('Optimum force with chosen magnitude and chosen piezo position from previous iteration')
    xlabel('Element')
    
 %}   

end

