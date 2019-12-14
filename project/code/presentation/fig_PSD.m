% clear
% close all


load('figures/PSD_64QAM_M1024_m2.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 

%fbmc
load('figures/PSD_FBMC_M1024.mat')%, 'fFBMC', 'sumFBMCSpec')

% hFig1 = figure('Position', figposition([46 15 30 30]));
figure;
plot(fOFDM,pow2db(psdOFDM),'linewidth',2);
hold on
plot(fCos,pow2db(psdCos),'linewidth',2);
plot(fFBMC-0.5,10*log10(sumFBMCSpec),'linewidth',2);

hold off
grid on
% xlim([-0.5, 0.5])
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)')
ylabel('Power/frequency (dBW/Hz)')
l=legend('OFDM','CMFB','FBMC');
set(gca, 'fontsize',13)
set(l, 'fontsize', 12, 'location', 'south');
set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 6 4])



save2file=0;
if(save2file)
saveas(gcf, 'figures/PSD_64QAM_M1024_m2.png' )
saveas(gcf, 'figures/PSD_64QAM_M1024_m2.svg' )
saveas(gcf, 'figures/PSD_64QAM_M1024_m2.pdf' )
!pdfcrop figures/PSD_64QAM_M1024_m2.pdf figures/PSD_64QAM_M1024_m2.pdf
end
% title(['PSD Comparison (' num2str(numRBs*rbSize) ' Subcarriers)'])