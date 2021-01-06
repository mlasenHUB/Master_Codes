function [ D, msfFRF, msfCF, oed, oemb, oeb, rCF, fs  ] = suppeffectA3(s, idloc1, idloc2,f, f2, Ma,Mpp,Mnp,Mnn,Mpn, VMMs,PP, OPT,r, label,bo, Nb )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


Dus=Ma(:,s,1);%Un-Suppresses Displacement
DsFRF= VMMs(:,s);%Suppressed Displacement from the FRF with optimum force f2
fs=PP(:,s)-f';
aux=OPT(s,4);%index to know where to find the displacement
n=length(fs);

if aux==1 
    DsCF= Mpp(:,s,idloc1,idloc2);
end
if aux==2 
    DsCF= Mnp(:,s,idloc1,idloc2);
end
if aux==3 
    DsCF= Mnn(:,s,idloc1,idloc2);
end
if aux==4 
    DsCF= Mpn(:,s,idloc1,idloc2);
end

[msfFRF]=MSF(DsFRF,Dus,bo, Nb);
[msfCF]=MSF(DsCF,Dus,bo, Nb);

%find the location of the maximum Over Excited Displacement and in how many
%blades it is present
raux=r;
raux(raux==1)=nan; %neglect unchanged values
oed=max(abs(raux(:,s)));
if oed<1
    oed=0;%no overexcitation et all
end
maxoed=(oed-1)*100;
oemb=0; %number of overexcited to the max level blades(to the maximum level)
oeb=0;
for i=1:n
    
    if r(i,s)==oed
        oemb=oemb+1;
    end
    if r(i,s)==-oed
        oemb=oemb+1;
    end
    if r(i,s)>1 || r(i,s)<-1
        oeb=oeb+1;
    end
end




D=[Dus , DsFRF , DsCF];%Displacements 
rCF= DsCF./Dus; %ratio for the displacement with Chosen Force wrt the unsuppressed displacente
rCFm=max(rCF);
%{
figure('units','normalized','outerposition',[0 0 1 1])
suptitle(label)
subplot(3,1,1)
    bar(D)
    xlabel('DOF')
    ylabel('Displacement')
    legend('UNSUPPRESED Displacement','SUPPRESSED WITH OPTIMUM FORCE', 'SUPPRESSED WITH A CHOSEN FORCE')
    title(['Mode: ', num2str(s),'- MSF TYPE B: ',num2str(msfFRF),'; MSF TYPE A/Max local MSF ', num2str(msfCF),' / ',num2str(rCFm) '; Max Overexcitation: ', num2str(maxoed), '; Number of Overexcited Blades (Max/Total): ', num2str(oemb),' / ', num2str(oeb)])

subplot(3,1,2)
    stem(fs)
    legend('Chosen froce')
    xlabel('Element')
    
subplot(3,1,3)

    stem(f2)
    hold on
    stem(f)
    legend('Force to total suppression', 'Excitation froce')
    xlabel('Element')
    
%}

end

