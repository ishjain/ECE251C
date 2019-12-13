%% F-OFDM vs. OFDM Modulation
% clc
% clf
close all

% clearvars

s = rng(211);       % Set RNG state for repeatability

%% System Parameters
% Define system parameters for the example. These parameters can be
% modified to explore their impact on the system.

numFFT = 1024;           % Number of FFT points
switch numFFT
    case 1024
        numRBs = 50;             % Number of resource blocks
        rbSize = 12;             % Number of subcarriers per resource block
        cpLen = 72;              % Cyclic prefix length in samples
    case 16
        numRBs = 5;             % Number of resource blocks
        rbSize = 2;             % Number of subcarriers per resource block
        cpLen = 4;              % Cyclic prefix length in samples
end
bitsPerSubCarrier = 6;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
% snrdB = 18;              % SNR in dB

for iter=1:1
snrlist = 18;%1:2:20;
for si = 1:length(snrlist)
snrdB = snrlist(si);
clear hh ff USout fout hout RxSymbolsCos
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator( ...
    'ModulationOrder', 2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');




% Generate data symbols
numDataCarriers = numRBs*rbSize;    % number of data subcarriers in subband
bitsIn = randi([0 1], bitsPerSubCarrier*numDataCarriers, 1);
symbolsIn = qamMapper(bitsIn);

%% OFDM Tx
%         takes symbolsIn (600x1) and returns symbolsIn (1096x1)
% Pack data into an OFDM symbol
offset = (numFFT-numDataCarriers)/2; % for band center
symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
    zeros(numFFT-offset-numDataCarriers,1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));

% Prepend cyclic prefix
txSigOFDM = [ifftOut(end-cpLen+1:end); ifftOut];

%% Cosine Modulated TMUX Transmit Processing
% TMUX parameter definition
m = 2;
N = 2*m*numFFT-1;
% Prototype filter definition
load(sprintf('prototype (M = %d, m = %d).mat',numFFT,m))
p0 = p0/sqrt(2*numFFT); % The prototype filter
if iscolumn(p0)
    p0 = p0.';
end
% TMUX (direct form) definition
n = 0:N;
for idx = 1:numFFT
    k = idx-1;
    thetak = (-1)^k*pi/4;
    hh(idx,:) = 2*p0.*cos(pi/numFFT*(k+1/2)*(n-N/2)+thetak);
    ff(idx,:) = 2*p0.*cos(pi/numFFT*(k+1/2)*(n-N/2)-thetak);
end
% Transmit data generation
symbolsInCos = [symbolsIn;zeros(numFFT-numDataCarriers,1)];
for idx = 1:numFFT
    USout(idx,:) = upsample(symbolsInCos(idx,:),numFFT); % The outputs of the upsamplers of all channels
    fout(idx,:) = conv(USout(idx,:),ff(idx,:)); % The outputs of the synthesis filters of all channels
end
txSigCos = sum(fout);
txSigCos = txSigCos(1:numFFT*2*m); % Keep only the nonzero samples

%have three symbols to address leakage issue.
txSigCos = [txSigCos zeros(1,2*numFFT)] + ...
    [zeros(1,numFFT) txSigCos zeros(1,numFFT)] + ...
    [zeros(1,2*numFFT) txSigCos];
% txSigCos = txSigCos(1,1:4*numFFT);
txSigCos = txSigCos(1,numFFT+(1:4*numFFT));

%% Channel

% Add WGN
rxSig = awgn(txSigOFDM, snrdB, 'measured');
rxSigCos = awgn(txSigCos, snrdB, 'measured');

rxSigCos = [0,rxSigCos];
% rxSigCos = [0,txSigCos]; % Just a delay

%% OFDM Receiver
% takes rxSig (1096x1) and return dataRxSymbols (600x1)
% Remove cyclic prefix
rxSymbol = rxSig(cpLen+1:end);

% Perform FFT
RxSymbols = fftshift(fft(rxSymbol));

% Select data subcarriers
dataRxSymbols = RxSymbols(offset+(1:numDataCarriers));


%% Cosine Modulated TMUX Receiver
for idx = 1:numFFT
    hout(idx,:) = conv(rxSigCos,hh(idx,:));
    RxSymbolsCos(idx,:) = downsample(hout(idx,:),numFFT);
end
dataRxSymbolsCos = RxSymbolsCos(1:numDataCarriers,5); % Select the 5th column to account for the delay



% Plot received symbols constellation
switch bitsPerSubCarrier
    case 2  % QPSK
        refConst = qammod((0:3).', 4, 'UnitAveragePower', true);
    case 4  % 16QAM
        refConst = qammod((0:15).', 16,'UnitAveragePower', true);
    case 6  % 64QAM
        refConst = qammod((0:63).', 64,'UnitAveragePower', true);
    case 8  % 256QAM
        refConst = qammod((0:255).', 256,'UnitAveragePower', true);
end
% constDiagRx = comm.ConstellationDiagram( ...
%     'ShowReferenceConstellation', true, ...
%     'ReferenceConstellation', refConst, ...
%     'Position', figposition([20 15 30 40]), ...
%     'EnableMeasurements', true, ...
%     'MeasurementInterval', length(dataRxSymbols), ...
%     'Title', 'OFDM Demodulated Symbols', ...
%     'Name', 'OFDM Reception', ...
%     'XLimits', [-1.5 1.5], 'YLimits', [-1.5 1.5]);
% constDiagRx(dataRxSymbols);
% 
% constDiagRx2 = comm.ConstellationDiagram( ...
%     'ShowReferenceConstellation', true, ...
%     'ReferenceConstellation', refConst, ...
%     'Position', figposition([20 15 30 40]), ...
%     'EnableMeasurements', true, ...
%     'MeasurementInterval', length(dataRxSymbols), ...
%     'Title', 'Cos Demodulated Symbols', ...
%     'Name', 'Cos Reception', ...
%     'XLimits', [-1.5 1.5], 'YLimits', [-1.5 1.5]);
% constDiagRx2(dataRxSymbolsCos);


% Channel equalization is not necessary here as no channel is modeled

%% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxBits = qamDemod(dataRxSymbols);
ber = BER(bitsIn, rxBits);

%% BER for Cos
% Perform hard decision and measure errors
rxBitsCos = qamDemod(dataRxSymbolsCos);
berCos = BER(bitsIn, rxBitsCos);

%% Display and plots
disp(['OFDM Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB) ' dB']);
disp(['COS Reception, BER = ' num2str(berCos(1)) ' at SNR = ' ...
    num2str(snrdB) ' dB']);

beralliter(si,:,iter) = [ber(1), berCos(1)];
end
end
berall = mean(beralliter,3);
save2file=0;
if(save2file)
save('figures/ber-snr_64QAM_M16_m2.mat','berall', 'snrlist');
end
%%plot BER

figure;
semilogy(snrlist, berall(:,1))
hold on
semilogy(snrlist, berall(:,2))
l=legend('OFDM', 'Cos TMux');
xlabel('SNR dB'); ylabel('BER')

%% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprCosTMUX] = PAPR(txSigCos.');
disp(['Peak-to-Average-Power-Ratio for Cosine modulated TMUX = ' num2str(paprCosTMUX) ' dB']);

%% OFDM Modulation with Corresponding Parameters
%
% For comparison, we review the existing OFDM modulation technique, using
% the full occupied band, with the same length cyclic prefix.

% Compute peak-to-average-power ratio (PAPR)
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = PAPR2(txSigOFDM);
disp(['Peak-to-Average-Power-Ratio for OFDM = ' num2str(paprOFDM) ' dB']);


%% Compare the PSD of transmit signals of OFDM and Cosine Modulated TMUX
% Normalize the average power for plotting
txSigOFDMPlt = txSigOFDM/sqrt(mean(abs(txSigOFDM.^2)));
txSigCosPlt = txSigCos/sqrt(mean(abs(txSigCos.^2)));
% Plot power spectral density (PSD) for OFDM signal
[psdOFDM,fOFDM] = periodogram(txSigOFDMPlt, rectwin(length(txSigOFDM)), 20000, ...
    1, 'centered');
[psdCos,fCos] = periodogram(txSigCosPlt, rectwin(length(txSigCos)), ...
    20000, 1, 'centered');

if(save2file)
   save('figures/PSD_64QAM_M16_m1.mat', 'psdCos', 'psdOFDM', 'fOFDM', 'fCos') 
end

hFig1 = figure('Position', figposition([46 15 30 30]));
plot(fOFDM,pow2db(psdOFDM));
hold on
plot(fCos,pow2db(psdCos));
hold off
grid on
% xlim([-0.5, 0.5])
axis([-0.5 0.5 -120 10]);
xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)')
ylabel('Power/frequency')
legend('OFDM','Cosine modulated TMUX')
title(['PSD Comparison (' num2str(numRBs*rbSize) ' Subcarriers)'])

% Restore RNG state
rng(s);