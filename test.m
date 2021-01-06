a=f*V;
N=length(a)/2;
Mmin2=zeros(4, 2*N);
for k = 1:2*N%Mode Shape
    for ii =1:2*N%Forcing pos
        Mmin2(ii, k)= -a(k)/V(ii, k);
    end    
end