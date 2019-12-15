clear
close all


load('figures/ber-snr_64QAM_M1024_m2.mat') %,'berall', 'snrlist');
ber2=berall;
load('figures/ber-snr_64QAM_M16_m2.mat') %,'berall', 'snrlist');
ber1=berall;

figure;
semilogy(snrlist, ber1(:,1), '-d' , 'MarkerSize',8,'linewidth',2)
hold on
semilogy(snrlist, ber2(:,1),'-o' ,  'MarkerSize',8,'linewidth',2)
grid on;


semilogy(snrlist, ber1(:,2), '-d' , 'MarkerSize',8,'linewidth',2)
hold on
semilogy(snrlist, ber2(:,2),'-o' ,  'MarkerSize',8,'linewidth',2)
grid on;

l=legend('M=16 OFDM', 'M=1024 OFDM','M=16 Cos', 'M=1024 Cos');
xlabel('SNR (dB)'); ylabel('BER')
set(gca, 'fontsize',13)
set(l, 'fontsize', 12);
set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 6 4])


saveas(gcf, 'figures/ber-snr_64QAM_impactM_m2.png' )
saveas(gcf, 'figures/ber-snr_64QAM_impactM_m2.svg' )

%%
load('figures/PSD_64QAM_M16_m2.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 
x1=fCos;
y1 = psdCos;

load('figures/PSD_64QAM_M1024_m2.mat');%, 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 
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
l=legend('M=16', 'M=1024');
set(gca, 'fontsize',13)
set(l, 'fontsize', 12, 'location', 'south');
set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 6 4])


save2file=0;
if(save2file)
saveas(gcf, 'figures/PSD_64QAM_impactM_m2.png' )
saveas(gcf, 'figures/PSD_64QAM_impactM_m2.svg' )
saveas(gcf, 'figures/PSD_64QAM_impactM_m2.pdf' )
!pdfcrop figures/PSD_64QAM_impactM_m2.pdf figures/PSD_64QAM_impactM_m2.pdf
end
% title(['PSD Comparison (' num2str(numRBs*rbSize) ' Subcarriers)'])