%close all;


% Point properties
num_of_points = N;
row_of_points = 2;

P=zeros(row_of_points, num_of_points);
P(1,:)=Vnd(:,5)';
P(2,:)=Vnb(:,5)';
% Create generic labels
P_labels = cell(num_of_points, 1);

for ii = 1:num_of_points
    P_labels{ii} = sprintf('Blisk %i', ii);
end

% Figure properties
figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);

% Axes properties
axes_interval = 1;
axes_precision = 1;

% Spider plot
spider_plot(P, P_labels, axes_interval, axes_precision,...
    'Marker', 'o',...
    'LineStyle', '-',...
    'LineWidth', 2,...
    'MarkerSize', 5);

% Title properties
title('Sample Spider Plot',...
    'Fontweight', 'bold',...
    'FontSize', 12);

% Legend properties
legend('show', 'Location', 'southoutside');
