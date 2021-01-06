function [ r ] = ratio( M1, M2 )
%Plot the ratio of each value for different matrices
%  M1= numerator.

r=M1./M2;
rabs=abs(r);

%figure
colormap(flipud(gray))
imagesc(rabs)
colorbar
xlabel('Modes')
ylabel('Element')
end

