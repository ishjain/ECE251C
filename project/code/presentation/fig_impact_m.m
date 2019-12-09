

load('figures/PSD_64QAM_M16_m1.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 
x1=fCos;
y1 = psdCos;

load('figures/PSD_64QAM_M16_m2.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 
x2=fCos;
y2 = psdCos;

% hFig1 = figure('Position', figposition([46 15 30 30]));
figure;
plot(x1,pow2db(y1),'linewidth',2);
hold on
plot(x2,pow2db(y2),'linewidth',2);


hold off
grid on
% xlim([-0.5, 0.5])
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)')
ylabel('Power/frequency (dBW/Hz)')
l=legend('m=1', 'm=2');
set(gca, 'fontsize',13)
set(l, 'fontsize', 12, 'location', 'south');

% saveas(gcf, 'figures/PSD_64QAM_impactM_m2.png' )
% saveas(gcf, 'figures/PSD_64QAM_impactM_m2.svg' )