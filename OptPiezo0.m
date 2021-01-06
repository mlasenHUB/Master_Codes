function [ PP, MSFf, idloc ] = OptPiezo0( f, Ma, Mp, Mn,s, bo, Nb )
%Finds the location to place a suppression force equal to p%(CHOSEN FORCE) of the
%excitation for (in/out of phase), to suppress all the possible modes. The
%optimum criterion per mode is the minimum MSF

%   s= mode to be suppressed
%   f2= optimum suppressing force(to suppress a mode to 0)
%   idloc= location of the optimum DOF to place the oftimum force in f2
%   VMMs= displacement obtained from the FRF suppresses with f2% 
%   bo=blades only?(1 yes, 0 no)


N=length(f);
mauxp=zeros(N,N);
mauxn=zeros(N,N);
for i=1:N
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

for i=1:N
   vaux=MP(:,i);
   vaux(vaux>min(vaux))=0;
   MP(:,i)=vaux;
end
for i=1:N
   vaux=MN(:,i);
   vaux(vaux>min(vaux))=0;
   MN(:,i)=vaux;
end

[rowp, colp]=find(MP>0);
[rown, coln]=find(MN>0);


 PP=zeros(N, N);%Piezo position optimum per mode and the corresponding force
 MSFf=zeros(N, N);%minimum obtainable MSF and where, per mode
for i=1:N
   if MP(rowp(i),i)<MN(rown(i),i)
       PP(rowp(i),i)= 0.01;
       MSFf(rowp(i),i)= MP(rowp(i),i);
   else
       PP(rown(i),i)= -0.01;
       MSFf(rown(i),i)= MN(rown(i),i);
   end
end

%find the optimum for mode s


[idloc]=find(MSFf(:,s)~=0); 






end

