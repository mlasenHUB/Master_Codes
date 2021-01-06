function MAC=mac(Phi1,Phi2, N, Nb)
% This function calculates mac between phi1 and phi2
% plot mac matrix
n=N+Nb;
MAC=zeros(n, n);
for j=1:n
for i = 1:n
    MAC(i,j)=Mac(Phi1(:,i),Phi2(:,j));
end
end

function mAc=Mac(Phi1,Phi2)
% This function calculates mac between phi1 and phi2
mAc= (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));
end


%figure
%bar3(MAC)
%title('MAC')

%figure
colormap(flipud(gray))
imagesc(MAC)
colorbar
title('MAC')
set(gca,'YDir','normal')

end