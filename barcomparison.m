function [ ] = barcomparison( c, Vsup, Vnsup, label1, label2 )
%   Plot the suppression effect in bars
% y3 displays the merged Vs and Vns.
%first bar in plot represent Vsup vectors and 2nd represents Vnsup
figure('units','normalized','outerposition',[0 0 1 1])
%suptitle(['Suppression effect for each mode and element' ])
n=length(Vsup(:,1));
bb=n/c;


y3=[];

for i=1:n
y2=[];%clean the vector of comparison y2 every iteration
for j=1:n
    y=[Vsup(j,i) Vnsup(j,i) ];
    y2= [y2 ; y];
    
    subplot(bb,c,i)       
    bar(y2)
    label=(['At freq. of mode: ', num2str(i)]);
    title(label)
    xlabel('Element')
    ylabel('Displacement')
end
legend(label1, label2)
y3=[y3, y2];

end  

%{
for i=1:n

subplot(bb,c,i)       
bar(y2)
legend('Suppressed', 'Unsuppressed')
label=(['At freq. of mode: ', num2str(i)]);
title(label)
xlabel('Element')
ylabel('Displacement')
end
%}

end

