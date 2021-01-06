function [ PP, MSFf, idloc ] = OptPiezoA2( Ma, Mp, Mn,Fp, Fn, s, N,bo,Nb )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n=N+Nb;
l=n+N;
mauxp=zeros(l,n);
mauxn=zeros(l,n);
for i=1:l
    [ms]=MSF( Mp(:,:,i),Ma(:,:,i), bo, Nb);
    [msn]=MSF( Mn(:,:,i),Ma(:,:,i),bo, Nb);
    mauxp(i,:)=ms;
    mauxn(i,:)=msn;
end

%minmsp=min(mauxp);
%minmsn=min(mauxn);

MP=mauxp;
MN=mauxn;
%MN(e, :)=nan; %trivial case of counteracting the excitation

for i=1:n
   vaux=MP(:,i);
   vaux(vaux>min(vaux))=0;
   a=nonzeros(vaux);b=length(a);
        if b>1
            for k=1:l
                if vaux(k)==a(1) && b>1 %eliminate duplicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end
   MP(:,i)=vaux;%Just the minimum MSF
end
for i=1:n
   vaux=MN(:,i);
   vaux(vaux>min(vaux))=0;
   a=nonzeros(vaux);b=length(a);
        if b>1
            for k=1:l
                if vaux(k)==a(1) && b>1 %eliminate duplicates in case there is more than one minimum value
                    vaux(k)=0;
                    b=b-1;
                end
            end
        end
   MN(:,i)=vaux;%Just the minimum MSF
end
  

[rowp, colp]=find(MP>0);%The only non-zero value is the minimum MSF.
[rown, coln]=find(MN>0);


 PP=zeros(n,n);%Piezo position optimum per mode and the corresponding force
 MSFf=zeros(l, n);%minimum obtainable MSF and where, per mode
 vf=0.01;
for i=1:n
   if MP(rowp(i),i)<= MN(rown(i),i)
       MSFf(rowp(i),i)= MP(rowp(i),i);
       PP(:,i)=Fp(rowp(i),:);
   end
   if MN(rown(i),i)<= MP(rowp(i),i)
       MSFf(rown(i),i)= MN(rown(i),i);
       PP(:,i)=Fn(rowp(i),:);
   end
end

%find the optimum for mode s
[idloc]=find(MSFf(:,s)~=0);

idloc=idloc(1);%in case there are more than 1 symmetric positions.


end

