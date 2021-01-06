function [ Ma, Mpp, Mnp, Mnn,Mpn, Fpp, Fnp,Fnn, Fpn, f2A ] = SuppressgivenforceA3( V,nr,p,d,f,N,Nb,td)
%Finds the suppression effect on a mode with a given force for Approach
%3(A3)
%   V= normalised eigen vectors
%   p= percentage of the excitacion force
%   d= squared eigen values
%   f= excitation force
%   Ma= Displacement without suppression (at resonances)
%   Mp= Displacement with supp force IN-PHASE with excitation
%   Mn= Displacement with supp force OUT-PHASE with excitation
%   td= tip to disk piezos (1= operating, 0=not in operation)


a=f*V;
n=N+Nb;
l=n+N;

M=zeros(n,n,l,l);
Mpp=zeros(n,n,l,l);
Mnp=zeros(n,n,l,l);
Mnn=zeros(n,n,l,l);
Mpn=zeros(n,n,l,l);

f2A=zeros(l,n);
%matrixes with all the possible combination of forces
Fpp=zeros(l,n,l);%vector forces in rows
Fnp=zeros(l,n,l);
Fnn=zeros(l,n,l);
Fpn=zeros(l,n,l);
%af=a+f2*V;%Piezo forcing in a certain frequency
for k = 1:n %frequency
    
    den=1i*nr(k)*d(k);%complex denominator at resonance
    for jj =1:l %%location of the suppression piezo 1(out of 2)
        for kk=1:l %% location of the piezo #2
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
            f2A(jj,:)=f2a;%forces in rows
                       
            fppp=f+f2A(jj,:)+f2A(kk,:);
            fpnp=f+f2A(jj,:)-f2A(kk,:);
            fpnn=f-(f2A(jj,:)-f2A(kk,:));
            fppn=f-(f2A(jj,:)+f2A(kk,:));
            
            Fpp(jj,:,kk)=fppp;
            Fnp(jj,:,kk)=fpnp;
            Fnn(jj,:,kk)=fpnn;
            Fpn(jj,:,kk)=fppn;
            
            
            app=fppp*V;
            anp=fpnp*V;
            ann=fpnn*V;
            apn=fppn*V;
            
            for ii=1:n  %DOF's displacement

                for m=1:n %mode contribution
                    M(ii,k,jj,kk)=M(ii,k,jj,kk)+ (V(ii,m)*a(m))/(d(m)-d(k)+den);   %wothout suppression force
                    Mpp(ii,k,jj,kk)= Mpp(ii,k,jj,kk)+ (V(ii,m)*app(m))/(d(m)-d(k)+den);%k: mode, ii:force position.(a(k):'pseudo' modal constant
                    Mnp(ii,k,jj,kk)= Mnp(ii,k,jj,kk)+ (V(ii,m)*anp(m))/(d(m)-d(k)+den);
                    Mnn(ii,k,jj,kk)= Mnn(ii,k,jj,kk)+ (V(ii,m)*ann(m))/(d(m)-d(k)+den);
                    Mpn(ii,k,jj,kk)= Mpn(ii,k,jj,kk)+ (V(ii,m)*apn(m))/(d(m)-d(k)+den);
                end
            end
        end
    end   
end
Ma=abs(M).*sign(sin(deg2rad(angle(M))));
Mpp=abs(Mpp).*sign(sin(deg2rad(angle(Mpp))));
Mnp=abs(Mnp).*sign(sin(deg2rad(angle(Mnp))));
Mnn=abs(Mnn).*sign(sin(deg2rad(angle(Mnn))));
Mpn=abs(Mpn).*sign(sin(deg2rad(angle(Mpn))));



end

