function [ Ma, Mp, Mn ] = Suppressgivenforce( V,nr,p,d,f)
%Finds the suppression effect on a mode whith a given force
%   V= normalised eigen vectors
%   p= percentage of the excitacion force
%   d= squared eigen values
%   f= excitation force
%   Ma= Displacement without suppression (at resonances)
%   Mp= Displacement with supp force IN-PHASE with excitation
%   Mn= Displacement with supp force OUT-PHASE with excitation

a=f*V;
N=length(a);
M=zeros(N,N,N);
Mp=zeros(N,N,N);
Mn=zeros(N,N,N);
%af=a+f2*V;%Piezo forcing in a certain frequency
for k = 1:N %frequency
    %f2a=zeros(1, 2*N);
    den=1i*nr(k)*d(k);%complex denominator at resonance
    for jj =1:N %%location of the suppression piezo
        for ii=1:N  %DOF
            f2a=zeros(1,N);
            
            f2a(jj)=(p/100)*max(abs(f));% a percentage of the max nonzero value in the excitation force
            a2=f2a*V;
            a3=a+a2;
            a4=a-a2;
            for m=1:N %mode contribution
            M(ii,k,jj)=M(ii,k,jj)+ (V(ii,m)*a(m))/(d(m)-d(k)+den);   %wothout suppression force     
            Mp(ii,k,jj)= Mp(ii,k,jj)+ (V(ii,m)*a3(m))/(d(m)-d(k)+den);%k: mode, ii:force position.(a(k):'pseudo' modal constant
            Mn(ii,k,jj)= Mn(ii,k,jj)+ (V(ii,m)*a4(m))/(d(m)-d(k)+den);
            end
        end
    end   
end
Ma=abs(M).*sign(sin(deg2rad(angle(M))));
Mp=abs(Mp).*sign(sin(deg2rad(angle(Mp))));
Mn=abs(Mn).*sign(sin(deg2rad(angle(Mn))));

end

