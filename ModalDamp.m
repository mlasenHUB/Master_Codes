function [ Hha, Hhta, ph ] = ModalDamp( b, g, r, V,d,  f, f2,w, ra, s,fdom )
%UNTITLED3 Summary of this function goes here
%Response with modal damping
%Hha= absolute values of mode contributing.
%Hhta= absolute total transfer function.
%b=proportional coefficient- Stiffness
%g=proportional coefficient- Mass
%r=Response point
%V=Normalized eigen vectors.
%d=eigen values
%f=Excitation Force.
%f2=Piezo Force.
%w=Window or not.
%ra= frequency range to apply the window.
%fdom=Frequency domain.


%Transfer function Multi forcing
a=f*V;
%x=NUM/DEN

N=length(a);
l=length(fdom);

if f2==0
    f2=zeros(1,N);
end
af=a+f2*V;%Piezo forcing in a certain frequency
if f==0
    f=zeros(1,N);
end

for i=1:N
    
    NUM(i)=V(r,i)*a(i);%Normalized eigenvectors V
end
%Modal Constant with forcing
for i=1:N
    
    NUMf(i)=V(r,i)*af(i);%Normalized eigenvectors V %These SHOULD NOT BE NORMALIZE FOR 
end
%C=bK+gM
nr=b+g./d;%damping loss factor d IS NOT the vector of damped natural frequencies.

for j=1:N
    for i=1:length(fdom)
        denhd=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        Hh(i,j)=NUM(j)/denhd;%Mode contribution SPD:Simple Prop. Damp.
        if w==1 && fdom(i)>=(d(s)^0.5-ra) && fdom(i)<=(d(s)^0.5+ra)
        Hh(i,j)=NUMf(j)/denhd;%Mode contribution SPD:Simple Prop. Damp.
        end
        if w==0 
        Hh(i,j)=NUMf(j)/denhd;%Mode contribution SPD:Simple Prop. Damp.
        end
    end
end
Hha=abs(Hh);%Absolute part

for i=1:length(fdom)
    Hhtc(i)=sum((Hh(i,:)));%Complex total transfer function.
end

for i=1:length(fdom)
    Hhta(i)=abs(Hhtc(i));%Absolute part of the total transfer function.
end

for i=1:length(fdom)
    ph(i)=angle(Hhtc(i));%Phase part of the total complex transfer function.
end

%ph=unwrap(ph);
ph=ph*180/pi;

end

