function [ msfvT,msfvM, minmsfT,minmsfM , RT, RM,  f2mT, f2mM, f2MinT, f2MinM, VTTs, VTTns, VMMs, VMMns, PHTs, PHTns,PHMs, PHMns ,y3T,y3M , smsfT, smsfM] = RunMCsim( nf, N, c,dT, dM , fdom, VT, VM, f, w, ra, s, b, g, plotnf, bo )
% Run Monte Carlo Simulation for several iterations
%  nf = number of random forces tried
% plot=1 then plot results
l=2*N;
msfvTa=zeros(1,nf);%vector
msfvT=zeros(1,nf);%vector
msfvMa=zeros(1,nf);%vector
msfvM=zeros(1,nf);%vector
f2mT=zeros(l,nf);%matrix
f2mM=zeros(l,nf);%matrix
tic
for i=1:nf
    i    
    [ msfTa, msfT, msfMa, msfM, f2T, f2M ] = MC( s,bo);
    
    msfvTa(i)=msfTa;
    msfvT(i)=msfT;
    msfvMa(i)=msfMa;
    msfvM(i)=msfM;
    f2mT(:,i)= f2T;
    f2mM(:,i)= f2M;
end



if bo==1 %blades only in the msf
    
    %force
    %tuned
    [minmsfT, ind]=min(abs(msfvT));
    f2MinT=f2mT(:,ind);
    %mistuned
    [minmsfM, indM]=min(abs(msfvM));
    f2MinM=f2mM(:,ind);
    
    %statistics
    %tuned
    xmT= mean(msfvT);
    stT= std(msfvT);
    vaT= var(msfvT);
    maxnrT=max(msfvT);   
    
    [minnrT, indT]=min(abs(msfvT));    
    statsT = [ xmT stT vaT minnrT maxnrT ];
    
    %mistuned
    xmM= mean(msfvM);
    stM= std(msfvM);
    vaM= var(msfvM);
    maxnrM=max(msfvM);   
    
    [minnrM, indTM]=min(abs(msfvM));    
    statsM = [ xmM stM vaM minnrM maxnrM ];
    
    %load msf from the optimal single DOF forces
    %tuned
    file=(['A1m', num2str(s)]); 
    smsfT=load(file,'msfT');
    smsfT=struct2cell(smsfT);
    smsfT=cell2mat(smsfT);
    %mistuned
    smsfM=load(file,'msfM');
    smsfM=struct2cell(smsfM);
    smsfM=cell2mat(smsfM);
    
    
else
    %force
    %tuned
    [minmsfT, ind]=min(abs(msfvTa));
    f2MinT=f2mT(:,ind);
    %mistuned
    [minmsfM, ind]=min(abs(msfvMa));
    f2MinM=f2mM(:,ind);
    
    %statistics
    %tuned
    xmT= mean(msfvT);
    stT= std(msfvT);
    vaT= var(msfvT);
    maxnrT=max(msfvT);   
    
    [minnrT, indT]=min(abs(msfvT));    
    statsT = [ xmT stT vaT minnrT maxnrT ];
    
    %mistuned
    xmM= mean(msfvM);
    stM= std(msfvM);
    vaM= var(msfvM);
    maxnrM=max(msfvM);   
    
    [minnrM, indTM]=min(abs(msfvM));    
    statsM = [ xmM stM vaM minnrM maxnrM ];
    
    %load msf from the optimal single DOF forces
    %tuned
    file=(['A1m', num2str(s)]); 
    smsfT=load(file,'msfTa');
    smsfT=struct2cell(smsfT);
    smsfT=cell2mat(smsfT);
    %mistuned
    smsfM=load(file,'msfMa');
    smsfM=struct2cell(smsfM);
    smsfM=cell2mat(smsfM);
    
    
end

%Plot the results for the min found from MC study for the tuned and
%mistuned sstem

%% Suppression effect on the tuned model
[ VTs, VTns, PHTs, PHTns, abslambda, xrTs, xtrTs, phrTs, xrTns, xtrTns, phrTns ] = SuppComparison( 1,dT, fdom, VT, f, f2MinT', w, ra, s, b, g );%Damping set to 1
[ VTTs, VTTns, y2T, y3T ] = suppeffect(VTs, VTns, 1, PHTs, PHTns );

RT=VTTs./VTTns;

%% Suppression effect on the mistuned model
[ VMs, VMns, PHMs, PHMns, abslambda, xrMs, xtrMs, phrMs, xrMns, xtrMns, phrMns ] = SuppComparison( 1,dM, fdom, VM, f, f2MinM', w, ra, s, b, g );%Damping set to 1
[ VMMs, VMMns, y2M, y3M ] = suppeffect(VMs, VMns, 1, PHMs, PHMns );
RM=VTTs./VTTns;

%% 
if plotnf==1
PlotTFpoints( c, xrTs, xtrTs, phrTs, fdom, 'TUNED-SUPPRESSED')
PlotTFpoints( c, xrTns, xtrTns, phrTns, fdom, 'TUNED-UNSUPPRESSED')
barcomparison( c, VTTs, VTTns, 'Tuned-Suppressed', 'Tuned-Unsuppressed');

PlotTFpoints( c, xrMs, xtrMs, phrMs, fdom, 'MISTUNED-SUPPRESSED')
PlotTFpoints( c, xrMns, xtrMns, phrMns, fdom, 'MISTUNED-UNSUPPRESSED')
barcomparison( c, VMMs, VMMns, 'Mistuned-Suppressed', 'Mistuned-Unsuppressed');

figure 
title(['Min MSF-T: ', num2str(minmsfT), '; Min MSF-M: ', num2str(minmsfM)])
stem(f, 'k-o');
xlim([1 2*N])
xlabel('DOF')
hold on
stem(f2MinT, 'r-x')
stem(f2MinM, 'b-x')
legend('Excitation Force','Optimal Tuned Suppressing Force', 'Optimal Mistuned Suppressing Force')
end

%% Statistics
figure('units','normalized','outerposition',[0 0 1 1])

%tuned
h1= subplot(2,1,1);
nbins=3*nf; %3 std devs right and left
%histfit(minmsfv,nbins)
histogram(msfvT,nbins)
grid on
title(['Monte Carlo simulation- Tuned System -Targeted mode: ', num2str(s)])
xlabel(['Msf for ' , num2str(nf)])
ylabel('Distribution of forces randomly generated')
h = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',h)
axes(h1)
line([smsfT,smsfT],[0,1],get(h1,'YLim'),'Color','red','LineStyle','--','LineWidth',2,'DisplayName','MSF with optimal single DOF force')
legend('Histogram of MSF', 'MSF with single DOF optimal force')

%mistuned
h2=subplot(2,1,2);
nbins=3*nf; %3 std devs right and left
%histfit(minmsfv,nbins)
histogram(msfvM,nbins)
grid on
title(['Monte Carlo simulation- Mistuned System -Targeted mode: ', num2str(s)])
xlabel(['Min msf for ' , num2str(nf)])
ylabel('Distribution of forces randomly generated')
h = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',h)
axes(h2)
line([smsfM,smsfM],[0,1],get(h2,'YLim'),'Color','red','LineStyle','--','LineWidth',2,'DisplayName','MSF with optimal single DOF force')
legend('Histogram of MSF', 'MSF with single DOF optimal force')

end

