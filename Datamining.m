%% DATA MINING: HERE I EXTRACT THE MOST IMPORTANT FEATURES FOR THE RESULTS OF THE THESIS.
datamining=6;%1:threshold study, 2: colateral effect, 3:opt for 1ppm bar plot, 4:ppm=1 only msf 1%
%% GET THE MISTUNING THRESHOLD AND PRINT IT
if datamining==1
    modes=12;
    
    
    figure('units','normalized','outerposition',[0 0 1 1])
    changek=zeros(modes, modes);
    for s =1:modes %modes to search
        file=(['k2s', num2str(s)]);%k-mistuning
        fsTTstruct=load(file,'fsTT'); %opt force for dif level of mist (1-10%)
        fsTTcell=struct2cell(fsTTstruct);
        fsTT=cell2mat(fsTTcell);
        
        fsMTstruct=load(file,'fsMT'); %opt force for dif level of mist (1-10%)
        fsMTcell=struct2cell(fsMTstruct);
        fsMT=cell2mat(fsMTcell);
        
        %check for a change in the opt force:
        
        aux1=fsTT-fsMT;
        aux2=max(abs(aux1));
        count=1;
        
        for i=1:length(aux2)
            
            if aux2(i)~=0
                aux2(i)=1;
                changek(count,s)=i;%count the times there is a change in the iptimal force
                count=count+1;
            end
        end
        aux3=[0,aux2];
        subplot(modes,2, 2*(s-1)+1)
        plot(aux3)
        xlim([1 11])
        xticks(1:1:11)
        xticklabels(0:1:10);
        
    end
    
    changem=zeros(modes, modes);
    for s =1:modes %modes to search
        file=(['m2s', num2str(s)]);%m-mistuning
        fsTTstruct=load(file,'fsTT'); %opt force for dif level of mist (1-10%)
        fsTTcell=struct2cell(fsTTstruct);
        fsTT=cell2mat(fsTTcell);
        
        fsMTstruct=load(file,'fsMT'); %opt force for dif level of mist (1-10%)
        fsMTcell=struct2cell(fsMTstruct);
        fsMT=cell2mat(fsMTcell);
        
        %check for a change in the opt force:
        
        aux1=fsTT-fsMT;
        aux2=max(aux1);
        count=1;
        
        for i=1:length(aux2)
            
            if aux2(i)~=0
                aux2(i)=1;
                changem(count,s)=i;%count the times there is a change in the iptimal force
                count=count+1;
            end
        end
        aux3=[0,aux2];
        subplot(modes,2, 2*(s))
        plot(aux3)
        xlim([1 11])
        xticks(1:1:11)
        xticklabels(0:1:10);
    end
%%    
    %%%%
    %%%% print forces before and after change
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('k-mistuning')
    
    for s =1:modes %modes to search
        file=(['k2s', num2str(s)]);%k-mistuning
        fsTTstruct=load(file,'fsTT'); %opt force for dif level of mist (1-10%)
        fsTTcell=struct2cell(fsTTstruct);
        fsTT=cell2mat(fsTTcell);
        
        fsMTstruct=load(file,'fsMT'); %opt force for dif level of mist (1-10%)
        fsMTcell=struct2cell(fsMTstruct);
        fsMT=cell2mat(fsMTcell);
        
        subplot(modes/2, 4, 2*(s-1)+1)
        drawforce(fsTT(:,1))%perfectly tuned force
        set(gca,'Visible','off')
        subplot(modes/2, 4, 2*(s))
        ind=changek(1,s);
        if ind~=0
            drawforce(fsMT(:,changek(1, s)))
            
        else
        end
        set(gca,'Visible','off')
    end
    
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('m-mistuning')
    for s =1:modes %modes to search
        file=(['m2s', num2str(s)]);%k-mistuning
        fsTTstruct=load(file,'fsTT'); %opt force for dif level of mist (1-10%)
        fsTTcell=struct2cell(fsTTstruct);
        fsTT=cell2mat(fsTTcell);
        
        fsMTstruct=load(file,'fsMT'); %opt force for dif level of mist (1-10%)
        fsMTcell=struct2cell(fsMTstruct);
        fsMT=cell2mat(fsMTcell);
        
        subplot(modes/2, 4, 2*(s-1)+1)
        drawforce(fsTT(:,1))%perfectly tuned force
        set(gca,'Visible','off')
        
        subplot(modes/2, 4, 2*(s))
        
        ind=changem(1,s);
        if ind~=0
            drawforce(fsMT(:,changem(1, s)))
            
        else
        end
        set(gca,'Visible','off')
    end
end
    
%% get the effect of one optimum in all other modes

if datamining==2
    %mistuned
    load('ColsOR12ats12M')
    
    x=1:1:12;
    
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('MISTUNED')
    for i=1:12%fixed optimal
        subplot(6,2,i)
        plot(msfcol(i,:))
        
        hold on
        plot(msfcurr(i,:))
        
        xlim([1 12])
        xticks(1:1:12)
        xticklabels(1:1:12);
        ylabel('MSF')
        legend(['Fixed optimal for mode: ',num2str(i)],'Local optimal','Location', 'southwest')
    end
  
    %tuned
    load('ColsOR12ats12T')
        
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('TUNED')
    for i=1:12
        subplot(6,2,i)
        plot(msfcol(i,:))
        hold on
        plot(msfcurr(i,:))
        xlim([1 12])
        xticks(1:1:12)
        xticklabels(1:1:12);
        ylabel('MSF')
        legend(['Fixed optimal for mode: ',num2str(i)],'Local optimal', 'Location', 'southwest')
    end
    
    %%%%
    %%%%
    
    %plot maxs ans mins
    load('ColsOR12ats12M')
    maxcol=zeros(1,12);
    mincol=zeros(1, 12);
    maxcur=zeros(1, 12);
    mincur=zeros(1, 12);
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('MISTUNED SYTEM')
    for i=1:12%fixed optimal
        for j=1:12%collateral modes 
            macol=max(bodcol(:,i,j));
            micol=min(bodcol(:,i,j));
            macur=max(bodcurr(:,i,j));
            micur=min(bodcurr(:,i,j));
            maxcol(j)=macol;
            mincol(j)=micol;
            maxcur(j)=macur;
            mincur(j)=micur;
        end
        ecolpos=maxcol-msfcol(i,:);
        ecolneg=mincol-msfcol(i,:);
        ecurpos=maxcur-msfcurr(i,:);
        ecurneg=mincur-msfcurr(i,:);
        subplot(6,2,i)
        plot(maxcol,'b-^')
        hold on
        plot(mincol,'b-v')
        plot(maxcur,'r-^')
        plot(mincur,'r-v')     
        
        xlim([1 12])
        xticks(1:1:12)
        xticklabels(1:1:12);
        ylabel('MSF')
    end
        legend('Max blade MSF, fixed piezo', 'Min blade MSF, fixed piezo','Max blade MSF, local optimal','Min blade MSF, local optimal')

  
    %tuned
    load('ColsOR12ats12T')
    maxcol=zeros(1,12);
    mincol=zeros(1, 12);
    maxcur=zeros(1, 12);
    mincur=zeros(1, 12);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle('TUNED SYSTEM')
    for i=1:12%fixed optimal
        for j=1:12%collateral modes 
            macol=max(bodcol(:,i,j));
            micol=min(bodcol(:,i,j));
            macur=max(bodcurr(:,i,j));
            micur=min(bodcurr(:,i,j));
            maxcol(j)=macol;
            mincol(j)=micol;
            maxcur(j)=macur;
            mincur(j)=micur;
        end
        ecolpos=maxcol-msfcol(i,:);
        ecolneg=mincol-msfcol(i,:);
        ecurpos=maxcur-msfcurr(i,:);
        ecurneg=mincur-msfcurr(i,:);
        subplot(6,2,i)
        plot(maxcol,'b-^')
        hold on
        plot(mincol,'b-v')
        plot(maxcur,'r-^')
        plot(mincur,'r-v')     
        
        xlim([1 12])
        xticks(1:1:12)
        xticklabels(1:1:12);
        ylabel('MSF')
    end
        legend('Max blade MSF, fixed piezo', 'Min blade MSF, fixed piezo','Max blade MSF, local optimal','Min blade MSF, local optimal')

  
    
end

%%
if datamining==3
    file=(['OPTs', num2str(12),'-ppm=1real'])';
    
    gy3T=load(file, 'gy3T');
    gy3T=struct2cell(gy3T);
    gy3T=cell2mat(gy3T);
    
    gy3M=load(file, 'gy3M');
    gy3M=struct2cell(gy3M);
    gy3M=cell2mat(gy3M);
    
    bb=3;
    c=4;
    modes=12;
    %tuned
    figure('units','normalized','outerposition',[0 0 1 1])
    
    
    for i=1:modes
        
        
        aux=gy3T(:, (2*(i-1)+1):(2*i));
        subplot(bb,c,i)
        bar(aux)
        %label=(['MSF:  ', num2str(msfT(i))]);
        %title(label)
        xlabel('Element')
        ylabel('Displacement')
    end
    suptitle('Tuned')
    
    %mistuned
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:modes
        
        
        aux=gy3M(:, (2*(i-1)+1):(2*i));
        subplot(bb,c,i)
        bar(aux)
        %label=(['MSF:  ', num2str(msfM(i))]);
        %title(label)
        xlabel('Element')
        ylabel('Displacement')
    end
    suptitle('Mistuned')
    
    
end
%%
if datamining==4
    msfTt=zeros(1,12);
    msfMt=zeros(1,12);
    
    for i=1:12
    file=(['OPTs', num2str(i),'-ppm=1real']);
    
    msfT=load(file, 'msfFRFT');
    msfT=struct2cell(msfT);
    msfT=cell2mat(msfT);
    
    msfTt(i)=msfT;
    
    msfM=load(file, 'msfFRFM');
    msfM=struct2cell(msfM);
    msfM=cell2mat(msfM);
    
    msfMt(i)=msfM;
    end
        
end


%%
if datamining==5
    idT=zeros(12,1);
    idM=zeros(12,1);
   for i=1:12 
   file=(['OPTs', num2str(i),'-ppm=1real']);
   
   idt=load(file, 'idlocTA2');
   idt=struct2cell(idt);
   idt=cell2mat(idt);
   
   idm=load(file, 'idlocMA2');
   idm=struct2cell(idm);
   idm=cell2mat(idm);
   
   idT(i)=idt;
   idM(i)=idm;
   end
end


%% datamine the tendency in the position of both piezos

if datamining==6
    
    
    s=1; %information of all the optimal is calculated for 1 mode only.    
        file=(['m3s', num2str(s)]);%k-mistuning or m-mistuning
        OPT=load(file,'OPT'); %positions of both piezos and msf for every mistuning level
        OPT=struct2cell(OPT);
        OPT=cell2mat(OPT);
        
        OPM=load(file,'OPM'); %positions of both piezos and msf for every mistuning level
        OPM=struct2cell(OPM);
        OPM=cell2mat(OPM);
        paT=zeros(12,10);
        pbT=zeros(12,10);
        paM=zeros(12,10);
        pbM=zeros(12,10);
        for mi=1:10
            for s=1:12 %modes, percentage of mistuning is inside evere archive
                
                paT(s,mi)=OPT(s,1,mi);
                pbT(s,mi)=OPT(s,2,mi);
                paM(s,mi)=OPM(s,1,mi);
                pbM(s,mi)=OPM(s,2,mi);
                
            end
        end
        
        %plot them for groups 1F-2F-3F
        figure('units','normalized','outerposition',[0 0 1 1])
        
        for i=9:12
        subplot(2,2, i-8)
        plot( paT(i,:),'b-')
        hold on
        plot(pbT(i,:),'b-')
        plot(paM(i,:),'rd')
        plot(pbM(i,:),'rd')
        lim=max([paT(i,:), pbT(i,:),paM(i,:),pbM(i,:)]);
        yticks(1:1:lim)
        ylim([1 lim]) 
        grid on
        xlabel('Mistuning [%]')
        ylabel('Optimal position')
        title(['Mode: #', num2str(i)])
        end
        legend('Tuned','Tuned','Mistuned','Mistuned')
    
end


%% change location of the excitation to check for a change in the suppression patterns

if datamining==7
    
    
    file=('m4sall');% m-mistuning
    PPT=load(file,'PPT'); %optimal forces
    PPT=struct2cell(PPT);
    PPT=cell2mat(PPT);
    
    PPM=load(file,'PPM'); %optimal forces
    PPM=struct2cell(PPM);
    PPM=cell2mat(PPM);
    
    f=load(file,'f');
    f=struct2cell(f);
    f=cell2mat(f);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    for s=1:12
        subplot(6, 4, 2*(s-1)+1)
        
        drawforce(PPT(:,s)-f')%perfectly tuned force
        set(gca,'Visible','off')
        
        subplot(6, 4, 2*(s))
        
        drawforce(PPM(:,s)-f')
        
        set(gca,'Visible','off')
    end
    
    
end

%% check the tendency of the piezos positions for a tuned system with a totally symmetric force
% in the tips

if datamining==8
    
    
    file=('Ts1');% m-mistuning
    PPT=load(file,'PPT'); %optimal forces
    PPT=struct2cell(PPT);
    PPT=cell2mat(PPT);
    
    PPM=load(file,'PPM'); %optimal forces
    PPM=struct2cell(PPM);
    PPM=cell2mat(PPM);
    
    f=load(file,'f');
    f=struct2cell(f);
    f=cell2mat(f);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    for s=1:12
        subplot(6, 2, s)
        
        drawforce(PPT(:,s)-f')%perfectly tuned force
        set(gca,'Visible','off')
     end
    
    
end






