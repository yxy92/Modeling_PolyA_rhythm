%phasescatterTotal.m
%scatter plot for phases



%% load data
clear;
load("BatchRunResultTable.mat");



%% extract PAR, relative amplitude greater than 0.2
PAR = result(result.Atrsc>=0.2 & result.AdeA>=0.2 & result.Adgrd>=0.2,:);
%% randomly sample 10% of the orignal data set
[m,n] = size(PAR);
y = randsample(m,1000);

PAR = PAR(y,:);

%% plot
figure;
colormap parula

% color scheme 
color1 = [73 0 146]/255;  % Trsc
color2 = [0 109 219]/255; % Dgrd
color3 = [146 73 0]/255;  % DeA



set(gcf,'DefaultLineLineWidth', 2)


ax=subplot(2,2,1);
ax.Position = [0.16 .68 0.28 0.28];
X1=scatter(PAR.Ptrsc,PAR.PL,10,color1,'filled');
pbaspect([1 1 1])
xlim([0 24])
ylim([0 24]);
xticks([0 12 24]);
yticks([0 12 24]);
set(gca,'LineWidth',2,'TickLength',[0.05 0.05]);
ax.FontSize = 20;

ax=subplot(2,2,2);
ax.Position = [0.6 0.68 0.28 0.28];
X4=scatter(PAR.Pdgrd,PAR.PL,10,color2,'filled');
pbaspect([1 1 1])
xlim([0 24])
ylim([0 24]);
xticks([0 12 24]);
yticks([0 12 24]);
set(gca,'LineWidth',2,'TickLength',[0.05 0.05]);
ax.FontSize = 20;

ax=subplot(2,2,3);
ax.Position = [0.16 0.15 0.28 0.28];
X2=scatter(PAR.PdeA,PAR.PL,10,color3,'filled');
pbaspect([1 1 1])
xlim([0 24])
ylim([0 24]);
xticks([0 12 24]);
yticks([0 12 24]);
set(gca,'LineWidth',2,'TickLength',[0.05 0.05]);
ax.FontSize = 20;



%save this figure
saveas(gcf,"PLPhaseScatter.png");