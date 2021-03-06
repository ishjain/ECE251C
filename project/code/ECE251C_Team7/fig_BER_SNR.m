clear
close all


load('figures/ber-snr_64QAM_M1024_m2.mat') %,'berall', 'snrlist');

%fbmc
load('figures/ber-snr_FBMC_64QAM_M1024.mat')%,'berFBMC', 'snrlistFBMC');

figure;
semilogy(snrlist, berall(:,1), '-d' , 'MarkerSize',8,'linewidth',2)
hold on
semilogy(snrlist, berall(:,2),'-o' ,  'MarkerSize',8,'linewidth',2)
hold on
semilogy(snrlistFBMC, berFBMC,'-h' ,  'MarkerSize',8,'linewidth',2)
grid on;


l=legend('OFDM', 'CMFB', 'FBMC');
xlabel('SNR (dB)'); ylabel('BER')
set(gca, 'fontsize',13)
set(l, 'fontsize', 12);
set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 6 4])


save2file=0;
if(save2file)
saveas(gcf, 'figures/ber-snr_64QAM_M1024_m2.png' )
saveas(gcf, 'figures/ber-snr_64QAM_M1024_m2.svg' )
saveas(gcf, 'figures/ber-snr_64QAM_M1024_m2.pdf' )
!pdfcrop figures/ber-snr_64QAM_M1024_m2.pdf figures/ber-snr_64QAM_M1024_m2.pdf
end
