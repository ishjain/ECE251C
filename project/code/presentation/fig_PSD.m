clear
close all


load('figures/PSD_64QAM_M1024_m2.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 


% hFig1 = figure('Position', figposition([46 15 30 30]));
figure;
plot(fOFDM,pow2db(psdOFDM));
hold on
plot(fCos,pow2db(psdCos));
hold off
grid on
% xlim([-0.5, 0.5])
axis([-0.5 0.5 -120 10]);
xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)')
ylabel('Power/frequency (dBW/Hz)')
l=legend('OFDM','Cos TMUX');
set(gca, 'fontsize',13)
set(l, 'fontsize', 12, 'location', 'south');

saveas(gcf, 'figures/PSD_64QAM_M1024_m2.png' )
saveas(gcf, 'figures/PSD_64QAM_M1024_m2.svg' )
% title(['PSD Comparison (' num2str(numRBs*rbSize) ' Subcarriers)'])