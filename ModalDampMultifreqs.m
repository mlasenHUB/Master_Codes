function [ Hha, Hhta, ph ] = ModalDampMultifreqs( b, g, r, V,d,  f, f2,w, ra, s,fdom, c )
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
% c = suppression frequency excitation proportional constant
% fs=c*fdom, where,(f2=f2(cos(fs*t))


%Transfer function Multi forcing
a=f*V;
a2=f2*V;
%x=NUM/DEN

N=length(a);
l=length(fdom);

if f2==0
    f2=zeros(1,N);
end

if f==0
    f=zeros(1,N);
end

for i=1:N
    
    NUMe(i)=V(r,i)*a(i);%Numerator of the excitation force
    NUMs(i)=V(r,i)*a2(i);%Numerator of the suppression force
end
%Modal Constant with forcing
%C=bK+gM
nr=b+g./d;%damping loss factor d IS NOT the vector of damped natural frequencies.

for j=1:N
    for i=1:length(fdom)
        denhde=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        denhd2=complex((d(j)-(c*fdom(i))^2), nr(j)*d(j));%denominator.
        
        Hh(i,j)=NUMe(j)/denhde;%Mode contribution SPD:Simple Prop. Damp.
        Hhs(i,j)=NUMs(j)/denhd2;%Mode contribution SPD:Simple Prop. Damp.
        
        %if w==1 && fdom(i)>=(d(s)^0.5-ra) && fdom(i)<=(d(s)^0.5+ra)
        %Hh(i,j)=NUMf(j)/denhde;%Mode contribution SPD:Simple Prop. Damp.
        %end
        if w==0 
        Hht(i,j)=Hh(i, j)+Hhs(i, j);%Mode contribution SPD:Simple Prop. Damp.
        end
    end
end
Hha=abs(Hht);%Absolute part

for i=1:length(fdom)
    Hhtc(i)=sum((Hht(i,:)));%Complex total transfer function.
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

