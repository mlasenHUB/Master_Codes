function [ Vs, Vns, PHs, PHns, abslambda, xr, xtr, phr, xrn, xtrn, phrn ] = SuppComparison( DAMP, d, fdom, V, f, f2, w, ra, s, b, g )
%% DISPLACEMENT OF MODES WITH PIEZO SUPPRESSION with and w/o Damping

%Find the response for all the possible positions

l=length(fdom);
n=length(d);
xr=cell(1, n);
xrn=cell(1, n);
xtr=zeros(l, n);
xtrn=zeros(l, n);
phr=zeros(l, n);
phrn=zeros(l, n);
%WITH SUPPRESSION (f2=Piezo Force =!0)
if DAMP==0
    
for i=1:n
    [ x, xt, ph ] = NoDamp( i, V,d,  f, f2, w, ra,s, fdom );
    xr{i}=x;
    xtr(:,i)=xt;  %WITH SUPPRESSION  (f2=!0) 
    phr(:,i)=ph;
end
%WITHOUT SUPPRESSION (f2=Piezo Force =0)
for i=1:n
    [ x, xt, ph ] = NoDamp( i, V,d,  f, 0, w, ra,s, fdom );
    xrn{i}=x;
    xtrn(:,i)=xt;  %WITHOUT SUPPRESSION (f2=0)
    phrn(:,i)=ph;
    
end

end

if DAMP==1
    
for i=1:n
    [ x, xt, ph ] = ModalDamp( b,g, i, V,d,  f, f2,w, ra, s, fdom );
    xr{i}=x;
    xtr(:,i)=xt;  %WITH SUPPRESSION  (f2=!0) 
    phr(:,i)=ph;
end
    
%WITHOUT SUPPRESSION
for i=1:n
    [ x, xt, ph ] = ModalDamp( b,g, i, V,d,  f, 0,w, ra, s, fdom );
    xrn{i}=x;
    xtrn(:,i)=xt;  %WITHOUT SUPPRESSION  (f2=0)  
    phrn(:,i)=ph;
end
    
end
% Create and display the modes with and without the suppression.
for i=1:n

    if DAMP==0
        [val, indf]=min(abs((d(i)^0.5-fdom)));% finding the closest position in fdom to the natural frequencies
    end
    if DAMP==1
        %Imaginary:
        nr=b+g./d;
        lambda2=d.*(complex(1,nr));
        lambda=sqrt(lambda2);%Complex number
        abslambda=abs(lambda);
        [val, indf]=min(abs((abslambda(i)-fdom)));
        indfv(i)=indf;
        %real
        %sd=sqrt(d);
        %chi=b.*sd/2+ g./(2.*sd);
        %dd=d.*(1-chi.^2); %Damped squared natural frequencies (not imaginary)
        %[val, indf]=min(abs((dd(i)^0.5-fdom)));
        
    end
    

vs=xtr(indf,:); %vs contains all the displacements of the nodes at one resonance frequency WITH SUPRESSION.
vns=xtrn(indf,:); %vns contains all the displacements of the nodes at one resonance frequency WITH NO SUPRESSION.

Vs(:,i)=vs;%Vs contains all the displacementes of the nodes at all resonant frequencies WITH SUPRESSION
Vns(:,i)=vns;%Vns contains all the displacementes of the nodes at all resonant frequencies WITH NO SUPRESSION

phs=phr(indf,:); % contains all the phases of the nodes at one resonance frequency WITH SUPRESSION.
phns=phrn(indf,:); % contains all the phases of the nodes at one resonance frequency WITH NO SUPRESSION.

PHs(:,i)=phs;%contains all the phasess of the nodes at all resonant frequencies WITH SUPRESSION
PHns(:,i)=phns;% contains all the phases of the nodes at all resonant frequencies WITH NO SUPRESSION
end

end

