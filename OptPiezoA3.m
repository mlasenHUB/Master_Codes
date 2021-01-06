function [ PP, MSFf, idloc1, idloc2, OPT, Dsup ] = OptPiezoA3( Ma, Mpp,Mnp,Mnn,Mpn,Fpp,Fnp,Fnn, Fpn,s, N,bo,Nb)
%FINDS THE OPTIMUM LOCATION OF 2 PIEZOS ACTING IN 2 DOF(EACH)
%FOR A GIVEN FORCE CONSIDERING THE MINIMUM MSF CRITERION
%   Detailed explanation goes here
n=N+Nb;
l=n+N;

mauxpp=zeros(l,n,l);
mauxnp=zeros(l,n,l);
mauxnn=zeros(l,n,l);
mauxpn=zeros(l,n,l);


for i=1:l%first piezo
    for j=1:l%second piezo,columns=MSF per mode
        [mspp]=MSF( Mpp(:,:,i,j),Ma(:,:,1,1),bo, Nb);
        [msnp]=MSF( Mnp(:,:,i,j),Ma(:,:,1,1),bo, Nb);
        [msnn]=MSF( Mnn(:,:,i,j),Ma(:,:,1,1),bo, Nb);
        [mspn]=MSF( Mpn(:,:,i,j),Ma(:,:,1,1),bo, Nb);
        
        mauxpp(i,:,j)=mspp;%msfs in columns for each mode
        mauxnp(i,:,j)=msnp;
        mauxnn(i,:,j)=msnn;
        mauxpn(i,:,j)=mspn;
        
    end
end

MPP=mauxpp;
MNP=mauxnp;
MNN=mauxnn;
MPN=mauxpn;

for i=1:n% frequency modes
    for j=1:l% second piezo(1st piezo in rows)
        vaux=MPP(:,i,j);
        vaux(vaux>min(vaux))=0;
        a=nonzeros(vaux);b=length(a);
        
        if b>1
            for k=1:l
                if vaux(k)==a(1) && b>1 %eliminate dulicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end
        MPP(:,i,j)=vaux;%Just the minimum MSF(there are many possible position to get the same MSF
    end
end
for i=1:n% frequency modes
    for j=1:l% second piezo(1st piezo in rows)
        vaux=MNP(:,i,j);
        vaux(vaux>min(vaux))=0;
        a=nonzeros(vaux);b=length(a);
        if b>1
            for k=1:l %each column 
                if vaux(k)==a(1) && b>1 %eliminate dulicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end
        MNP(:,i,j)=vaux;%Just the minimum MSF(there are many possible position to get the same MSF
    end
end
for i=1:n% frequency modes
    for j=1:l% second piezo(1st piezo in rows)
        vaux=MNN(:,i,j);
        vaux(vaux>min(vaux))=0;
        a=nonzeros(vaux);b=length(a);
        if b>1
            for k=1:l
                if vaux(k)==a(1) && b>1 %eliminate dulicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end
        MNN(:,i,j)=vaux;%Just the minimum MSF(there are many possible position to get the same MSF
    end
end
for i=1:n% frequency modes
    for j=1:l% second piezo(1st piezo in rows)
        vaux=MPN(:,i,j);
        vaux(vaux>min(vaux))=0;
        a=nonzeros(vaux);b=length(a);
        if b>1
            for k=1:l
                if vaux(k)==a(1) && b>1 %eliminate dulicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end  
        
        MPN(:,i,j)=vaux;%Just the minimum MSF(there are many possible position to get the same MSF)
    end
end

rowPP=zeros(l,n);
rowNP=zeros(l,n);
rowNN=zeros(l,n);
rowPN=zeros(l,n);




for j=1:l%second piezo
    %optimum position for the first piezo is in the number in the rows,
    %columns are the freqs of the modes and j is the second piezo
    [rowpp, col]=find(MPP(:,:,j)>0);%The only non-zero value is the minimum MSF.
    [rownp, col]=find(MNP(:,:,j)>0);
    [rownn, col]=find(MNN(:,:,j)>0);
    [rowpn, col]=find(MPN(:,:,j)>0);
    
        
    rowPP(j,:)=rowpp;%second piezo in rows number, first in cell, modes in columns
    rowNP(j,:)=rownp;
    rowNN(j,:)=rownn;
    rowPN(j,:)=rowpn;
    
end




PP=zeros(n,n);%Piezo position optimum per mode and the corresponding force
MSFf=zeros(l,l,n);%minimum obtainable MSF and where, per mode
vf=0.01;
for i=1:n % modes
    for j=1:l % second piezo
        
        if MPP(rowPP(j,i),i,j) <= MNP(rowNP(j,i),i,j) && MPP(rowPP(j,i),i,j)<= MNN(rowNN(j,i),i,j) && MPP(rowPP(j,i),i,j)<= MPN(rowPN(j,i),i,j)
            
            ind1=rowPP(j,i);%index of the optimum first piezo
            ind2=j;%index of the optimum second piezo
            MSFf(ind1,ind2,i)= MPP(rowPP(j,i),i,j);
            %obtain the absolute minimum(not just per position of the
            %second piezo)
        end
        if MNP(rowNP(j,i),i,j) <= MNN(rowNN(j,i),i,j) && MNP(rowNP(j,i),i,j)<= MPN(rowPN(j,i),i,j) && MNP(rowNP(j,i),i,j)<= MPP(rowPP(j,i),i,j)
            
            ind1=rowNP(j,i);%index of the optimum first piezo
            ind2=j;%index of the optimum second piezo
            MSFf(ind1,ind2,i)= MNP(rowNP(j,i),i,j);
        end
        if MNN(rowNN(j,i),i,j) <= MPN(rowPN(j,i),i,j) && MNN(rowNN(j,i),i,j)<= MPP(rowPP(j,i),i,j) && MNN(rowNN(j,i),i,j)<= MNP(rowNP(j,i),i,j)
            
            ind1=rowNN(j,i);%index of the optimum first piezo
            ind2=j;%index of the optimum second piezo
            MSFf(ind1,ind2,i)= MNN(rowNN(j,i),i,j);
        end
        if MPN(rowPN(j,i),i,j) <= MPP(rowPP(j,i),i,j) && MPN(rowPN(j,i),i,j)<= MNP(rowNP(j,i),i,j) && MPN(rowPN(j,i),i,j)<= MNN(rowNN(j,i),i,j) 
            
            ind1=rowPN(j,i);%index of the optimum first piezo
            ind2=j;%index of the optimum second piezo
            MSFf(ind1,ind2,i)= MPN(rowPN(j,i),i,j);
        end
    end
end


%Eliminate the redundant positions: ind1=ind2(piezos superposed)
for i=1:n
    pospm=MSFf(:,:,i);%positions per mode
    pospm(logical(eye(l)))=nan;
    MSFf(:,:,i)=pospm;
end



%find the optimum locations for mode s, ind1=piezo1, ind2=piezo2
for k=1:n %modes
    [maxv,inds1]=max(MSFf(:,:,k));
    maxv(maxv==0)=nan;
    [minv,ind2]=min(maxv);
    ind1=inds1(ind2);
    MSFf(:,:,k)=0;
    MSFf(ind1, ind2, k)=minv;
end


%find from which matrix is the minimum extracted and build the suppression
%force, also save the suppressed displacement.
mins=MSFf(MSFf>0); %all the minimums per mode
OPT=zeros(n,4);%matrix with opt piezo #1,2 ans min MSF for all the modes
Dsup=zeros(n,n);
for i=1:n %modes
    
    [ind1pp, ind2pp]=find(MPP(:,i,:)==mins(i));%locations of piezo #1 and #2
    [ind1np, ind2np]=find(MNP(:,i,:)==mins(i));
    [ind1nn, ind2nn]=find(MNN(:,i,:)==mins(i));
    [ind1pn, ind2pn]=find(MPN(:,i,:)==mins(i));
    
    if ~isempty(ind1pp) %if there is a result from the 'find' function.
        %build the optimum force and record the optimiums locations/msf
        ind1=ind1pp(1);
        ind2=ind2pp(1);
        OPT(i,1)=ind1;
        OPT(i,2)=ind2;
        OPT(i,3)=MPP(ind1,i,ind2);
        OPT(i,4)=1;%1 for PP
        PP(:,i)=Fpp(ind1,:,ind2);%forces to the columns per mode
        Dsup(:,i)=Mpp(:,i,ind1, ind2);%displacement of each DOF with suppression
    end
        
    if ~isempty(ind1np) %if there is a result from the 'find' function.
        %build the optimum force:
        ind1=ind1np(1);
        ind2=ind2np(1);
        OPT(i,1)=ind1;
        OPT(i,2)=ind2;
        OPT(i,3)=MNP(ind1,i,ind2);
        OPT(i,4)=2;%2 for NP
        PP(:,i)=Fnp(ind1,:,ind2);%forces to the columns per mode
        Dsup(:,i)=Mnp(:,i,ind1, ind2);%displacement of each DOF with suppression
    end
    if ~isempty(ind1nn) %if there is a result from the 'find' function.
        %build the optimum force:
        ind1=ind1nn(1);
        ind2=ind2nn(1);
        OPT(i,1)=ind1;
        OPT(i,2)=ind2;
        OPT(i,3)=MNN(ind1,i,ind2);
        OPT(i,4)=3;%3 for NN
        PP(:,i)=Fnn(ind1,:,ind2);%forces to the columns per mode
        Dsup(:,i)=Mnn(:,i,ind1, ind2);%displacement of each DOF with suppression
    end
    if ~isempty(ind1pn) %if there is a result from the 'find' function.
            
        
        %build the optimum force:
        ind1=ind1pn(1);%peacking one in case there are more
        ind2=ind2pn(1);
        
        OPT(i,1)=ind1;
        OPT(i,2)=ind2;
        OPT(i,3)=MPN(ind1,i,ind2);
        OPT(i,4)=4;%4 for PN
        PP(:,i)=Fpn(ind1,:,ind2);%forces to the columns per mode
        Dsup(:,i)=Mpn(:,i,ind1, ind2);%displacement of each DOF with suppression
    end
        
end


idloc1=OPT(s,1);
idloc2=OPT(s,2);
end

