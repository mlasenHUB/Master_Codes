function [ nfd ] = natfreqND( d, N, Nb )
%Plot the natural frequencies against the ND from the eigen problem
%   d= squared natural frequencies
%   N= number of diks 
%   Nb=Number of blades

nf=sqrt(d);
nfd=uniquetol(nf, 0.000001*nf(1));% nf without double modes
l=length(nfd);
n=l/(Nb/N+1);
fb=nfd(1:n);%first Bending family
sb=nfd(n+1:2*n);%second Bending family

%figure
plot(0:n-1,fb, '-o')
hold on
plot(0:n-1, sb, 'r-x')
if Nb>N
tb= nfd(2*n+1:l);%third bending family.
plot(0:n-1, tb, 'b-x')
end
hold off
set(gca,'xtick',0:n-1)

xlabel('Nodal diameters')
ylabel('Natural Frequencies')
legend('1-B', '2-B', '3-B')


end

