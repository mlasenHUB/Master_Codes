function [MACe, MACf, rSUS, nfT, nfM ] = GeneralandMAC( N,Nb, f, f2T, f2M, VM, VT, VMMM, VTTT, K, Kor, Ma, Maor, dT, dM, VMMs, VMMns )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
km=diag(K); km=km(N+1:end); %stiffness of the Mistuned blades
kt=diag(Kor); kt=kt(N+1:end); %stiffness of the Tuned blades
mm=diag(Ma); mm=mm(N+1:end);% mass of the Mistuned blades
mt=diag(Maor); mt=mt(N+1:end);% mass of the Tuned blades
n=N+Nb;
 
if length(f)>2*N
   nb=Nb/N;%number of blades per disk
   kM=diag(K); kmr=kM(Nb+1:end); %stiffness of the Mistuned blades RIGHT
   kml=kM(N+1:Nb); kml=kml-kmr; %stiffness of the Mistuned blades LEFT
   kT=diag(Kor); ktr=kT(Nb+1:end); %stiffness of the Mistuned blades RIGHT
   ktl=kT(N+1:Nb); ktl=ktl-ktr; %stiffness of the Mistuned blades LEFT
   
   mm=diag(Ma); mml=mm(N+1:Nb);mmr=mm(Nb+1:end);% mass of the Mistuned blades LEFT and RIGHT
   mt=diag(Maor); mtl=mt(N+1:Nb);mtr=mt(Nb+1:end);% mass of the Tuned blades LEFT and RIGHT
   
   figure('units','normalized','outerposition',[0 0 1 1])

    subplot(3,3,1)
    stem(kml, 'r-o')
    hold on
    stem(kmr, 'b-o')
    stem(ktl, 'r-x')
    stem(ktr, 'b-x')
    hold off
    xlabel('Blades')
    set(gca, 'XTick', 1:N)
    title('STIFFNESS: Mistuned: o - Tuned: x, Left/Up: red - Right/Down: blue')
    
    subplot(3,3,2)
    stem(mml, 'r-o')
    hold on
    stem(mmr, 'b-o')
    stem(mtl, 'r-x')
    stem(mtr, 'b-x')
    hold off
    xlabel('Blades')
    set(gca, 'XTick', 1:N)
    title('MASS: Mistuned: o - Tuned: x, Left: red - Right: blue')

   
end   

if length(f)<=2*N
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,3,1)
    stem(km, '-o')
    hold on
    stem(kt, 'x')
    hold off
    xlabel('Blades')
    set(gca, 'XTick', 1:N)
    title('STIFFNESS: Mistuned: o - Tuned: x')

    subplot(3,3,2)
    stem(mm, '-o')
    hold on
    stem(mt, 'x')
    hold off
    xlabel('Blades')
    set(gca, 'XTick', 1:N)
    title('MASS: Mistuned: o - Tuned: x')
end

subplot(3,3,3)
stem(f, '-o')
xlim([1 N+Nb])
xlabel('DOF')
title('Excitation Force: o ; Piezo Force: x(tuned/mistuned:red/blue')
hold on
stem(f2T, 'r-x')
stem(f2M, 'b-x')
hold off

%ratio of the mistuned to the tuned blade  for the corresponding ith nat frequencies, extracted from the FRF
subplot(3,3,4)
[ rMT ] = ratio( VMMM, VTTT ); title('Displacement comparison Mistuned/Tuned')
%MAC of the mode shapes from the eigen problem and extracted from the FRF
%Eigen 
subplot(3,3,5)
MACe=mac(VM,VT, N, Nb);xlabel('Tuned Modes (eig)');ylabel('Mistuned Modes (eig)')
%FRF
subplot(3,3,6)
MACf=mac(VMMM,VTTT, N, Nb);xlabel('Tuned ODS (FRF)');ylabel('Mistuned ODS (FRF)')
subplot(3,3,7)
[ nfT ] = natfreqND( dT, N, Nb); title('Tuned');set(gca,'ytick',linspace(min(nfT),max(nfT),15));grid on;ylim([min(nfT) max(nfT)])
subplot(3,3,8)
[ nfM ] = natfreqND( dM, N, Nb ); title('Mistuned');set(gca,'ytick',linspace(min(nfM),max(nfM),15));grid on; ylim([min(nfM) max(nfM)])

%ratio of the mistuned suppressed to unsuppressed displacement  for the corresponding ith nat frequencies, extracted from the FRF
subplot(3,3,9)
[ rSUS ] = ratio( VMMs, VMMns ); title('Displacement comparison Mistuned:Suppressed/Unsuppressed');set(gca,'ytick',1:N+Nb);set(gca,'xtick',1:N+Nb)

end

