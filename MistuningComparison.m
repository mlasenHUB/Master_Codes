function [ VTT, VMM, PHT, PHM, abslambdaT, abslambdaM, xrT, xtrT, phrT, xrM, xtrM, phrM ] = MistuningComparison( dT, dM, fdom, VT, VM, f, w, ra, s, b, g )
%% displacement with and without mistuning

%Find the response for all the possible positions
%WITH SUPPRESSION (f2=Piezo Force =!0)
l=length(fdom);
n=length(dT);
xrM=cell(1, n);
xrT=cell(1, n);
xtrT=zeros(l, n);
xtrM=zeros(l, n);
phrT=zeros(l, n);
phrM=zeros(l, n);

%Tuned
for i=1:n
    [ x, xt, ph ] = ModalDamp( b,g, i, VT,dT,  f, 0,w, ra, s, fdom );
    xrT{i}=x;
    xtrT(:,i)=xt;  %WITHOUT SUPPRESSION  (f2=0)  
    phrT(:,i)=ph;
end
%Mistuned
for i=1:n
    [ x, xt, ph ] = ModalDamp( b,g, i, VM, dM,  f, 0,w, ra, s, fdom );
    xrM{i}=x;
    xtrM(:,i)=xt;  %WITHOUT SUPPRESSION  (f2=0)  
    phrM(:,i)=ph;
end


% Create and display the modes with and without the mistuning.
for i=1:n

%Imaginary:
%Tuned
nrT=b+g./dT;%Loss factor

lambda2T=dT.*(complex(1,nrT));
lambdaT=sqrt(lambda2T);%Complex number
abslambdaT=abs(lambdaT);
[val, indfT]=min(abs((abslambdaT(i)-fdom)));
%indfvT(i)=indfT; %Indexes of the resonances in the domain
%Mistuned
nrM=b+g./dM;%Loss factor

lambda2M=dM.*(complex(1,nrM));
lambdaM=sqrt(lambda2M);%Complex number
abslambdaM=abs(lambdaM);
[val, indfM]=min(abs((abslambdaM(i)-fdom)));
%indfvM(i)=indfM; %Indexes of the resonances in the domain


%real
%sd=sqrt(d);
%chi=b.*sd/2+ g./(2.*sd);
%dd=d.*(1-chi.^2); %Damped squared natural frequencies (not imaginary)
%[val, indf]=min(abs((dd(i)^0.5-fdom)));
        

    

vT=xtrT(indfT,:); %vT contains all the displacements of the nodes at one resonance frequency TUNED
vM=xtrM(indfM,:); %vns contains all the displacements of the nodes at one resonance frequency MISTUNED.

VTT(:,i)=vT;%Vs contains all the displacementes of the nodes at all resonant frequencies TUNED
VMM(:,i)=vM;%Vns contains all the displacementes of the nodes at all resonant frequencies MISTUNED

phT=phrT(indfT,:); %vs contains all the displacements of the nodes at one resonance frequency TUNED
phM=phrM(indfM,:); %vns contains all the displacements of the nodes at one resonance frequency MISTUNED

PHT(:,i)=phT;%Vs contains all the displacementes of the nodes at all resonant frequencies TUNED
PHM(:,i)=phM;%Vns contains all the displacementes of the nodes at all resonant frequencies MISTUNED
end

end

