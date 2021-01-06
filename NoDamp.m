function [ x, xt, ph ] = NoDamp( r, V,d,  f, f2, w, ra, s, fdom )
%UNTITLED3 Summary of this function goes here
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

if norm(f2)==0
   f2=zeros(1,N);
end

af=a+f2*V;%Piezo forcing in a certain frequency
if f==0
    f=zeros(1,N);
end

NUM=zeros(1, N);
NUMf=zeros(1, N);
for i=1:N
    
    NUM(i)=V(r,i)*a(i);%Normalized eigenvectors V
end
%Modal Constant with forcing

for i=1:N
    
    NUMf(i)=V(r,i)*af(i);%Normalized eigenvectors V
end
x=zeros(l,N);
for j=1:N
    for i=1:l
        DEN=d(j)-fdom(i)^2;%denominator.
        x(i,j)=NUM(j)/DEN;%Mode contribution.
        %Windowing: the piezo acts just in a range ra of the frequency.
        if w==1 && fdom(i)>=(d(s)^0.5-ra) && fdom(i)<=(d(s)^0.5+ra)% && j==s
            
           x(i,j)=NUMf(j)/DEN;%Modes contribution.
        end
        if w==0 
            
           x(i,j)=NUMf(j)/DEN;%Modes contribution.
        end
    end
end

xt=zeros(1, l);
 for i=1:l
    xt(i)=sum((x(i,:)));%Total transfer function without damping.
 end

 ph=zeros(1, l);
 for i=1:l
    ph(i)=angle(xt(i));%Phase part of the total complex transfer function.
 end
 
%ph=unwrap(ph);
ph=ph*180/pi;

end

