function [ Ma, Mp, Mn, Fp, Fn, f2A ] = SuppressgivenforceA2( V,nr,p,d,f,N,Nb,td)

% Finds the suppression effect on a mode with a given force for Approach
% 2(A2)
%   V= normalised eigen vectors
%   p= percentage of the excitacion force
%   d= squared eigen values
%   f= excitation force
%   Ma= Displacement without suppression (at resonances)
%   Mp= Displacement with supp force IN-PHASE with excitation
%   Mn= Displacement with supp force OUT-PHASE with excitation

a=f*V;
n=N+Nb;
l=n+N;
M=zeros(n,n,l);
Mp=zeros(n,n,l);
Mn=zeros(n,n,l);

Fp=zeros(l,n); %forces in rows for l posssible positions (p=positive w.r.t. the excitation)
Fn=zeros(l,n); %forces in rows for l posssible positions (n=negative w.r.t the excitation)


f2A=zeros(l,n);%force


for k = 1:n %frequency
    den=1i*nr(k)*d(k);%complex denominator at resonance
    for jj =1:l %%location of the suppression piezo
        f2a=zeros(1,n);
        vf=(p/100)*max(abs(f));%(magnitude) a percentage of the nonzero value in the excitation force
        if jj<=2*N
            f2a(jj)=vf;
            f2a(jj+N)=-vf;
        end
        if jj>2*N && jj<=n
            f2a(jj)=vf*td;
            f2a(jj-2*N)=-vf*td;
        end
        if jj> n && jj<n+N
            f2a(jj-n)=vf;
            f2a(jj-n+1)=-vf;
        end
        if jj==n+N
            f2a(jj-n)=vf;
            f2a(1)=-vf;
        end
        
        f2A(jj,:)=f2a;
        
        fp=f+f2A(jj,:);
        fn=f-f2A(jj,:);
        
        Fp(jj,:)=fp;
        Fn(jj,:)=fn;
        
        a2=f2a*V;
        a3=a+a2;
        a4=a-a2;
        for ii=1:n  %DOF
            
            for m=1:n %mode contribution
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

