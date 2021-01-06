function [ SF, SFv ] = Suppression( R, s, N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Rs=R(:, s); % suppressed
Rss=Rs;
N=length(Rs);
Sup_disk = Rs(1:N/2);
Sup_blades = Rs(N/2+1:end);

d_min= min(Sup_disk);
d_max= max(Sup_disk);
b_min= min(Sup_blades);
b_max= max(Sup_blades);
r_d= d_min/d_max;
r_b=  b_min/b_max;

%Just considered the suppressed:

indices = find(abs(Rs)>1);
Rss(indices) = [];

m_all=mean(Rs);
m_s=mean(Rss); %mean of just the suppressed components


        

SF{1}= (['Most Suppressed Disk: ' , num2str(d_min)]);
SF{2}= (['Most Suppressed Blade: ', num2str(b_min)]);
SF{3}= (['Less Suppressed Disk: ', num2str(d_max)]);
SF{4}= (['Less Suppressed Blade: ', num2str(b_max)]);
SF{5}= (['Ratio of the most suppressed vs less suppressed of the disks: ', num2str(r_d)]);
SF{6}= (['Ratio of the most suppressed vs less suppressed of the blades: ', num2str(r_b)]);
SF{7}= (['Mean of the suppression of the components in the mode: ', num2str(m_all)]);
SF{8}= (['Mean of only the successfully suppressed components in the mode: ', num2str(m_s)]);

SFv=[d_min d_max b_min b_max r_d r_b m_all m_s];

end

