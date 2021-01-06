function [Mmin2, Mnom2] = Suppressfinddoubleforce( f, V, nbs,td )
%Find the optimal value of a force acting in TWO DOFs, in TENSION or
%COMPRESION, Valid for models with more than 1 DOF per Blade
%   nbs= number of blades in series per disk

a=f*V;
n=length(a);
N=n/(nbs+1);
Mmin2=zeros(n+N,n);

for k = 1:n%Mode Shape
    for ii =1:n+N %Forcing pos
        % piezo connected 'to the right': d(disk)->lb(left blade)->rb(right blade) in the model d+lb, lb+rb, rb+d
        if ii<=2*N
            Mmin2(ii, k)= -a(k)/(V(ii, k)-V(ii+N, k)); %d+lb, lb+rb,
        end
        if td
            if ii>2*N && ii<= n  % disk to tip connections
                Mmin2(ii, k)= -a(k)/(V(ii, k)-V(ii-2*N, k)); %rb+d
            end
        end
        if ii>n && ii<n+N %between disks connection
            Mmin2(ii,k)=-a(k)/(V(ii-n, k)-V(ii-n+1, k));
            
        end
        if ii==n+N
            Mmin2(ii,k)= -a(k)/(V(ii-n, k)-V(1, k)); %last disk connected with first disk
        end
    end    
end

Mmin2(~isfinite(Mmin2))=NaN;%Neglect infinite values of forcing
Mmin2(isnan(Mmin2))=0;
Mnom2=Mmin2./max(abs(Mmin2(:)));

end

