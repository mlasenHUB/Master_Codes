function [ f2 ] = CreateSuppForce(MT, MM , MISTK, MISTM, idloc, Pt , N,Nb,s, ppm)
%Create the suppressing force
%M= minimisation force matrix, for 1 or 2 piezos per mode to be suppressed
%MT=Tuned
%MM=Mistuned
%N= tblisk units
%Nb=number of blades per blisk unit
%s=mode you want to suppress
if ppm==1
    if MISTK==0 && MISTM==0 %if there is no mistuning use the tuned optimised matrix
        Pd= MT(idloc,s);
    end
    if MISTK ==1 || MISTM==1 %if there is mistuning use the mistuned optimised matrix
        Pd = MM(idloc, s);
    end
    
    f2d=zeros(1, Nb+N);
    
    pp1=idloc;
    
    %action force (as in action-reaction in the piezo)
    
    if pp1<= N+Nb
        f2d(pp1)=Pd;
    else
        f2d(pp1-(N+Nb))=Pd; %between disks piezos
    end
    
    %Create the reaction force depending on the postition
    if pp1<=2*N
        f2d(pp1+N)=-Pd;
    end
    if pp1>2*N && pp1<=N+Nb 
        f2d(pp1-2*N)=-Pd;
    end
    if pp1>N+Nb %between disk piezo
        f2d(pp1-(N+Nb)+1)=-Pd;
        if pp1==N+Nb+N%last disk's reaction in first disk
            f2d(1)=-Pd;
        end
    end
    
f2=f2d;
end




if ppm==2
    pp1=Pt(s,1);
    pp2=Pt(s,2);
    %both piezos have the same force
    if MISTK ==0 && MISTM==0 %if there is mistuning use the mistuned optimised matrix
        Pd = MT(pp1,pp2, s);
    end
    if MISTK ==1 || MISTM==1 %if there is mistuning use the mistuned optimised matrix
        Pd = MM(pp1,pp2, s);
    end
    %action force (as in action-reaction in the piezo)
    f2d=zeros(1, Nb+N);
    if pp1<= N+Nb
        
        f2d(pp1)=Pd;
    else
        f2d(pp1-(N+Nb))=Pd; %between disks piezos
    end
    if pp2<= N+Nb
        
        f2d(pp2)=f2d(pp2)+Pd;%not replace but add in case the superspose in one location
    else
        f2d(pp2-(N+Nb))=f2d(pp2-(N+Nb))+Pd; %between disks piezos
    end
    
    %Create the reaction force depending on the postition
    if pp1<=2*N
        f2d(pp1+N)=f2d(pp1+N)-Pd;
    end
    if pp1>2*N && pp1<=N+Nb 
        f2d(pp1-2*N)=f2d(pp1-2*N)-Pd;
    end
    if pp1>N+Nb %between disk piezo
        f2d(pp1-(N+Nb)+1)=f2d(pp1-(N+Nb)+1)-Pd;
        if pp1==N+Nb+N%last disk's reaction in first disk
            f2d(1)=f2s(1)-Pd;
        end
    end
    
    if pp2<=2*N
        f2d(pp2+N)=f2d(pp2+N)-Pd;
    end
    if pp2>2*N && pp2<=N+Nb 
        f2d(pp2-2*N)=f2d(pp2-2*N)-Pd;
    end
    if pp2>N+Nb %between disk piezo
        f2d(pp2-(N+Nb)+1)=f2d(pp2-(N+Nb)+1)-Pd;
        if pp2==N+Nb+N%last disk's reaction in first disk
            f2d(1)=f2d(1)-Pd;
        end
    end
       
    
f2=f2d;
    
end


end

