clear
close all


load('figures/ber-snr_64QAM_M1024_m2.mat') %,'berall', 'snrlist');

figure;
semilogy(snrlist, berall(:,1), '-bd' , 'MarkerSize',8,'linewidth',2)
hold on
semilogy(snrlist, berall(:,2),'-ro' ,  'MarkerSize',8,'linewidth',2)

grid on;


l=legend('OFDM', 'Cos TMux');
xlabel('SNR (dB)'); ylabel('BER')
set(gca, 'fontsize',13)
set(l, 'fontsize', 12);

saveas(gcf, 'figures/ber-snr_64QAM_M1024_m2.png' )
saveas(gcf, 'figures/ber-snr_64QAM_M1024_m2.svg' )

