
a=f*V;
F=1;
%Mmin

for k = 1:2*N
    f2=zeros(1, 2*N);
    for ii =1:2*N
        f2(ii)=F;
        a2=f2*V;
        M(k, ii)= a(k)+a2(ii);
    end
    
end

%Best positions to place the piezo to suppress that mode:
Mmin=zeros(4, 2*N);
for i =1:2*N
    [minval, pos]=min(abs(M(:,i)));
    Mmin(1, i)=pos;
    Mmin(2, i)=minval;
    Mmin(3, i)=abs(min(a(i)));
    if Mmin(2,i) < Mmin(3,i)
        Mmin(4, i)= 1;
    end
end

%Mmin2
%Best forcing value to place it at any position, to suppress that mode
Mmin2=zeros(4, 2*N);
for k = 1:2*N%Mode Shape
    for ii =1:2*N%Forcing pos
        Mmin2(ii, k)= -a(k)/V(ii, k);
    end    
end


