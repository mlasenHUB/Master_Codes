function [ msf ] = MSF( Vs, Vns, bo, Nb )
%UNTITLED Summary of this function goes here
%   Vs= suppressed matrix
%   Vns= unsuppressed matrix
%   bo= blades-only (1 if yes, 0 if no)
%   Nb=number of blades
[a,b]=size(Vs);
if bo==0
    
    if a==1 || b==1
        for i=1:max(a,b)
            numv(i)=Vs(i)*Vns(i);
            denv(i)=Vns(i)*Vns(i);
        end
        msf=sum(abs(numv))/sum(abs(denv));
    end
       
    if a>1 && b>1
        N=length (Vs(:,1));
                
        for j=1:N
            
            for i=1:N
                numv(i)=Vs(i,j)*Vns(i,j);
                denv(i)=Vns(i,j)*Vns(i,j);
            end
            
            msf(j)=sum(abs(numv))/sum(abs(denv));
        end
    end
end

if bo==1
    N=length (Vs(:,1));
    
    n=N-Nb;
    
    if a==1 || b==1
        for i=n+1:max(a,b)%just blades
            numv(i)=Vs(i)*Vns(i);
            denv(i)=Vns(i)*Vns(i);
        end
        msf=sum(abs(numv))/sum(abs(denv));
    end 
    
    if a>1 && b>1
        
        
        for j=1:N %modes
            
            for i=n+1:N %elements, just blades
                numv(i)=Vs(i,j)*Vns(i,j);
                denv(i)=Vns(i,j)*Vns(i,j);
            end
            
            msf(j)=sum(abs(numv))/sum(abs(denv));
        end
    end
    
end

end

