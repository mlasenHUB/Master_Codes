function [ Ft, Pt , Mf, idloc] = OptPiezo( Mmin3,s )
%Return the matrix of forces and positions to suppress with 1...N piezos
%   Detailed explanation goes here
Ma=abs(Mmin3);
[r,c]=size(Ma);%number of rows and columns (c=modes, r=positions)
%Local optimum tu suppres only mode s

Ma(Ma==0)=nan;%eliminate the disregarded cases of the DISK connected to the BASE-BLADE

[locmin,idloc]=min(Ma(:,s));


[m1, pm1]=min(sum(Ma,2));%minimum for 1 piezo
p11=Ma(pm1,:);

[min_val,idx]=min(Ma(:,:));

idx(idx==pm1)=nan;


Ft=zeros(r,c,r);%forces 
Ft(1,:,1)=p11;
Pt=zeros(r,r);
Pt(1,1)=pm1;
Mf=zeros(3*r,r);%matrix of maximum forces and amount of modes suppressed


for p=2:r
    
    
    Ft(:,:,p)=Ft(:, :,p-1);
    
    [pmp, frepp]=mode(idx);%position of the most repeated minimum for p groups
    
    Pt(p,:)= Pt(p-1,:);
    Pt(p,p)=pmp;
        
        for i=1:c
            if idx(i)== pmp
                Ft(1,i,p)=0;
                Ft(p,i,p)=Ma(pmp, i);
            end    
        end
    
    idx(idx==pmp)=nan; 
    
end


for i=1:r
    
    for j=1:i
        v=Ft(j,:,i);
        ma=max(v);%max force
        mi=min(v(v>0));%min force
        if isempty(mi)
            mi=0;
        end
        n=nnz(v);%numbers of modes suppress by the piezo
        Mf(3*i-2,j)=ma;
        Mf(3*i-1,j)=mi;
        Mf(3*i,j)=n;
    end
end

end

