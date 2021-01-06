function [ j,k ] = FOP( N, f,V, F )
%fop=Find Optimal Position
%f is the original excitation force.
%F is the value of the piezo force
%N=Degrees of Freedom(
%V is the matrix of modal shapes(in columns)
%   Finds F*V for every possible forcing position j and mode k, then it
%   calculates a factor M to minimise.

a=f*V;
for k = 1:2*N
    for ii =1:2*N
        if ii==k
            f2(ii)=F;
        else
            f2(ii)=0;
        end
        a2=f2*V;
        M(ii)= a(k)-a2(ii);
    end
    
end




end

