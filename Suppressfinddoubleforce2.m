function [Mmin2, Mnom2] = Suppressfinddoubleforce2( f, V, nbs, td )
%Find the optimal value of a force acting in TWO DOFs, in TENSION or
%COMPRESION, acting with  more than one piezo per mode.
%Valid for models with more than 1 DOF per Blade
%   nbs= number of blades in series per disk
%   td= tip-disk forces (0: not operating, 1: operating)

a=f*V;
n=length(a);

N=n/(nbs+1);
Mmin2=nan*ones(n+N,n+N,n);
for k = 1:n%Mode Shape
    % piezo connected 'to the right': d(disk)->lb(left blade)->rb(right blade) in the model d+lb, lb+rb, rb+d
    for ii =1:n+N%Forcing pos for the first piezo
        for jj=1:n+N%Forcing pos for the second piezo
            
            if ii<=2*N
                if jj<=2*N
                    Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii+N, k)+V(jj, k)-V(jj+N, k)); %d+lb, lb+rb,
                end
                if td
                    if jj>2*N && jj<= n
                        Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii+N, k)+V(jj, k)-V(jj-2*N, k));
                    end
                end
                if jj>n && jj<n+N %between disks connection
                    Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii+N, k)+V(jj-n, k)-V(jj-n+1, k));
                    
                end
                if jj==n+N
                    Mmin2(ii,jj,k)= -a(k)/(V(ii, k)-V(ii+N, k)+V(jj-n, k)-V(1, k)); %last disk connected with first disk
                end
            end
            if td
                if ii>2*N && ii<=n
                    if jj<=2*N
                        Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii-2*N, k)+V(jj, k)-V(jj+N, k)); %d+lb, lb+rb,
                    end
                    if jj>2*N && jj<= n
                        Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii-2*N, k)+V(jj, k)-V(jj-2*N, k));
                    end
                    if jj>n && jj<n+N %between disks connection
                        Mmin2(ii,jj, k)= -a(k)/(V(ii, k)-V(ii-2*N, k)+V(jj-n, k)-V(jj-n+1, k));
                    end
                    if jj==n+N
                        Mmin2(ii,jj,k)= -a(k)/(V(ii, k)-V(ii-2*N, k)+V(jj-n, k)-V(1, k)); %last disk connected with first disk
                    end
                end
            end
            
            if ii>n && ii<n+N %between disks piezo
                if jj<=2*N
                    Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(ii-n+1, k)+V(jj, k)-V(jj+N, k)); %d+lb, lb+rb,
                end
                if td
                    if jj>2*N && jj<= n
                        Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(ii-n+1, k)+V(jj, k)-V(jj-2*N, k));
                        
                    end
                end
                if jj>n && jj<n+N %between disks connection
                    Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(ii-n+1, k)+V(jj-n, k)-V(jj-n+1, k));
                end
                if jj==n+N
                    Mmin2(ii,jj,k)= -a(k)/(V(ii-n, k)-V(ii-n+1, k)+V(jj-n, k)-V(1, k)); %last disk connected with first disk
                end
            end
            if ii==n+N
                if jj<=2*N
                    Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(1, k)+V(jj, k)-V(jj+N, k)); %d+lb, lb+rb,
                    
                end
                if td
                    if jj>2*N && jj<= n
                        Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(1, k)+V(jj, k)-V(jj-2*N, k));
                        
                    end
                end
                if jj>n && jj<n+N %between disks connection
                    Mmin2(ii,jj, k)= -a(k)/(V(ii-n, k)-V(1, k)+V(jj-n, k)-V(jj-n+1, k));
                end
                if jj==n+N
                    Mmin2(ii,jj,k)= -a(k)/(V(ii-n, k)-V(1, k)+V(jj-n, k)-V(1, k)); %last disk connected with first disk
                end
            end
        end
    end
    
    Mmin2(~isfinite(Mmin2))=NaN;%Neglect infinite values of forcing
    %Mmin2(isnan(Mmin2))=0;
    Mnom2=Mmin2./max(abs(Mmin2(:)));
    
end

