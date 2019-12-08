% BERanalysisRx.m

%%--Call transmitter code
BERanalysisTx;

% Rx processing params
FFT_OFFSET                    = 0;          % Number of CP samples to use in FFT (on average)
SNR                           = 30;         %dB
%%--Future: Add a loop in SNR and plot BER-SNR curve

rx_vec_air_temp =  tx_vec_air ;


%% ------Add Channel-----
%  channel = zeros(1,1600);
channel(1) = 1;%*exp(1j*2*pi*rand);
%  channel(800) = 0.9*exp(1j*2*pi*rand);
%  channel(400) = 0.8*exp(1j*2*pi*rand);


rx_vec_air = conv(rx_vec_air_temp,channel);



%% --------Add noise---------

%%--Calculate noise power
sig_pow = mean(abs(rx_vec_air).^2);
SNR_num = db2pow(SNR);
noise_pow = sig_pow/SNR_num;
noise = 1/sqrt(2)*sqrt(noise_pow)*complex(randn(1,length(rx_vec_air)), randn(1,length(rx_vec_air)));

raw_rx_dec = rx_vec_air + noise;


%% Extract LTS for channel estimate

lts_ind=1;
payload_ind = 2.5*N_SC + lts_ind;

rx_lts = raw_rx_dec(lts_ind : lts_ind+2.5*N_SC-1);
rx_lts1 = rx_lts(-N_SC-FFT_OFFSET + [1.5*N_SC+1:2.5*N_SC]);
rx_lts2 = rx_lts(-FFT_OFFSET + [1.5*N_SC+1:2.5*N_SC]);

rx_lts1_f = fft(rx_lts1);
rx_lts2_f = fft(rx_lts2);

%%--Calculate channel estimate from average of 2 training symbols
rx_H_est = lts_f .* (rx_lts1_f + rx_lts2_f)/2;
rx_h_est = ifft(rx_H_est);
% figure; plot(abs(rx_h_est))

%if want to avoid channel estimate, just use
% rx_h_est = channel;
% rx_H_est = fft(rx_h_est);

%% Rx payload processing

% Extract the payload samples (integral number of OFDM symbols following preamble)
payload_vec = raw_rx_dec(payload_ind : end);

if(doFBMC) %doFBMC defined in BERanalysisTx.m
    
    % take payload_vec (e.g 1x(80x100) for OFDM) and return rx_syms (1x6400 for OFDM)
    for idx = 1:M
        hout(idx,:) = [0,conv(tx,h(idx,:))];
        xHat(idx,:) = downsample(hout(idx,:),M);
    end
    xRec = xHat(:,3:end-2);
    
    rx_syms=xRec;
    
    
else %do OFDM
    % Reshape
    payload_mat = reshape(payload_vec, (N_SC+CP_LEN), N_OFDM_SYMS);
    
    % Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
    payload_mat_noCP = payload_mat(CP_LEN-FFT_OFFSET+[1:N_SC], :);
    
    % Take the FFT
    syms_f_mat = fft(payload_mat_noCP, N_SC, 1);
    
    % Equalize (zero-forcing, just divide by complex chan estimates)
    syms_eq_mat = syms_f_mat ./ repmat(rx_H_est.', 1, N_OFDM_SYMS);
    
    % Extract useful symbols (removing pilot symbol and unloaded symbols)
    payload_syms_mat = syms_eq_mat(SC_IND_DATA, :);
    
    
    rx_syms = reshape(payload_syms_mat, 1, N_DATA_SYMS);
end


%% Demodulate
demod_fcn_bpsk = @(x) double(real(x)>0);
demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));
demod_fcn_64qam = @(x) (32*(real(x)>0)) + (16*(abs(real(x))<0.6172)) + (8*((abs(real(x))<(0.9258))&&((abs(real(x))>(0.3086))))) + (4*(imag(x)>0)) + (2*(abs(imag(x))<0.6172)) + (1*((abs(imag(x))<(0.9258))&&((abs(imag(x))>(0.3086)))));

switch(MOD_ORDER)
    case 2         % BPSK
        rx_data = arrayfun(demod_fcn_bpsk, rx_syms);
    case 4         % QPSK
        rx_data = arrayfun(demod_fcn_qpsk, rx_syms);
    case 16        % 16-QAM
        rx_data = arrayfun(demod_fcn_16qam, rx_syms);
    case 64        % 64-QAM
        rx_data = arrayfun(demod_fcn_64qam, rx_syms);
end

%% Plot Results
cf = 20;

% Rx signal
cf = cf + 1;
figure(cf); clf;
subplot(2,1,1);
plot(real(rx_vec_air), 'b');
% axis([0 length(rx_vec_air) -TX_SCALE TX_SCALE])
grid on;
title('Rx Waveform (I)');

subplot(2,1,2);
plot(imag(rx_vec_air), 'r');
% axis([0 length(rx_vec_air) -TX_SCALE TX_SCALE])
grid on;
title('Rx Waveform (Q)');



% Channel Estimates
cf = cf + 1;

rx_H_est_plot = repmat(complex(NaN,NaN),1,length(rx_H_est));
rx_H_est_plot(SC_IND_DATA) = rx_H_est(SC_IND_DATA);


x = (20/N_SC) * (-(N_SC/2):(N_SC/2 - 1));

figure(cf); clf;
subplot(2,1,1);
stairs(x - (20/(2*N_SC)), fftshift(real(rx_H_est_plot)), 'b', 'LineWidth', 2);
hold on
stairs(x - (20/(2*N_SC)), fftshift(imag(rx_H_est_plot)), 'r', 'LineWidth', 2);
hold off
axis([min(x) max(x) -1.1*max(abs(rx_H_est_plot)) 1.1*max(abs(rx_H_est_plot))])
grid on;
title('Channel Estimates (I and Q)')

subplot(2,1,2);
bh = bar(x, fftshift(abs(rx_H_est_plot)),1,'LineWidth', 1);
shading flat
set(bh,'FaceColor',[0 0 1])
axis([min(x) max(x) 0 1.1*max(abs(rx_H_est_plot))])
grid on;
title('Channel Estimates (Magnitude)')
xlabel('Baseband Frequency (MHz)')


%% Symbol constellation
cf = cf + 1;
figure(cf); clf;

plot(payload_syms_mat(:),'ro','MarkerSize',1);
axis square; axis(1.5*[-1 1 -1 1]);
grid on;
hold on;

plot(tx_syms_mat(:),'bo');
title('Tx and Rx Constellations')
legend('Rx','Tx','Location','EastOutside');




% EVM & SNR
cf = cf + 1;
figure(cf); clf;

evm_mat = abs(payload_syms_mat - tx_syms_mat).^2;
aevms = mean(evm_mat(:));
snr = 10*log10(1./aevms);

subplot(2,1,1)
plot(100*evm_mat(:),'o','MarkerSize',1)
axis tight
hold on
plot([1 length(evm_mat(:))], 100*[aevms, aevms],'r','LineWidth',4)
myAxis = axis;
h = text(round(.05*length(evm_mat(:))), 100*aevms+ .1*(myAxis(4)-myAxis(3)), sprintf('Effective SNR: %.1f dB', snr));
set(h,'Color',[1 0 0])
set(h,'FontWeight','bold')
set(h,'FontSize',10)
set(h,'EdgeColor',[1 0 0])
set(h,'BackgroundColor',[1 1 1])
hold off
xlabel('Data Symbol Index')
ylabel('EVM (%)');
legend('Per-Symbol EVM','Average EVM','Location','NorthWest');
title('EVM vs. Data Symbol Index')
grid on

subplot(2,1,2)
imagesc(1:N_OFDM_SYMS, (SC_IND_DATA - N_SC/2), 100*fftshift(evm_mat,1))

grid on
xlabel('OFDM Symbol Index')
ylabel('Subcarrier Index')
title('EVM vs. (Subcarrier & OFDM Symbol)')
h = colorbar;
set(get(h,'title'),'string','EVM (%)');
myAxis = caxis();
if (myAxis(2)-myAxis(1)) < 5
    caxis([myAxis(1), myAxis(1)+5])
end



%% Calculate Rx stats (SER, BER, EVM)

sym_errs = sum(tx_data ~= rx_data);
bit_errs = length(find(dec2bin(bitxor(tx_data, rx_data),8) == '1'));
rx_evm   = sqrt(sum((real(rx_syms) - real(tx_syms)).^2 + (imag(rx_syms) - imag(tx_syms)).^2)/(length(SC_IND_DATA) * N_OFDM_SYMS));

fprintf('\nResults:\n');
fprintf('Num Bytes:   %d\n', N_DATA_SYMS * log2(MOD_ORDER) / 8);
fprintf('Sym Errors:  %d (of %d total symbols)\n', sym_errs, N_DATA_SYMS);
fprintf('Bit Errors:  %d (of %d total bits)\n', bit_errs, N_DATA_SYMS * log2(MOD_ORDER));

%%--distance between figures
distFig;






