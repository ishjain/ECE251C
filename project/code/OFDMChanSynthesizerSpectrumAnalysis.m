% A channel synthesizers in various forms and the spectrum of the output
clc
clf
close all
clearvars

%% Channel Synthesizer parameter definition
M = 64; % The number of channels
L = 1e4; % The number of frames

%% Symbol definition
ch = 50; % The data-carrying channel
sigmar = 1;
sigmai = 1;
% rng(0);
Sr = sigmar*randn(1,L);
% rng(1);
Si = sigmai*randn(1,L);
S = zeros(M,L);
S(ch,:) = complex(Sr,Si);

%% Prototype filter definition
W = exp(-1j*2*pi/M);
f0 = ones(1,M)/sqrt(M); % The prototype filter

%% A channel synthesizer in direct form
% Loop over each channel
for idx = 1:M
    k = idx-1;
    f(idx,:) = f0.*W.^(-k*[0:M-1]); % The filter
    w(idx,:) = upsample(S(idx,:),M); % The symbols after upsampling
    p(idx,:) = conv(w(idx,:),f(idx,:)); % The filtered symbols
end
x = sum(p); % Sum over all channels

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 2 polyphase is used)
% % The polyphase components of the prototype filter go to each channel
% R = ones(M,1)/sqrt(M); 
% FF2 = dftmtx(M); % The DFT matrix
% 
% % % Compute DFT by looping over each frame
% % for idx = 1:L
% %     SS(:,idx) = FF2*S(:,idx); % Data symbols after DFT
% % end
% 
% SS = FF2*S;
% 
% % Shuffling
% SM = [zeros(M-1,1),eye(M-1);1,zeros(1,M-1)];
% SS = SM*SS;
% 
% % Loop over each channel for filtering
% for idx = 1:M
%     k = idx-1;
%     Q(idx,:) = conv(R(idx,:),SS(idx,:));
% end
% 
% % S/P conversion
% x2 = flipud(Q);
% x2 = x2(:).';

%% A channel synthesizer in polyphase form (version 2)
% % The polyphase components of the prototype filter go to each channel
% R = ones(M,1)/sqrt(M); 
% FF = conj(dftmtx(M)); % The IDFT matrix
% 
% % Compute IDFT by looping overSA.ShowLegend = true; each frame
% for idx = 1:L
%     SS(:,idx) = FF*S(:,idx); % Data symbols after IDFT
% end
% 
% % Loop over each channel for filtering
% for idx = 1:M
%     k = idx-1;
%     Q(idx,:) = conv(R(idx,:),SS(idx,:));
% end
% 
% % S/P conversion
% x3 = flipud(Q);
% x3 = x3(:).';

%% What about in the time domain
FF = conj(dftmtx(M))/sqrt(M); % The IDFT matrix
x4 = FF*S;
% Add CP
x4 = [x4(64-15:64,:);x4];
x4 = x4(:).';

%% Spectrum analyzer definition
fspace = 312.5e3; % Let's say we are using it as WiFi
fs = M*fspace;
fprintf('fs = %d MHz\n',fs/1e6)
SA = dsp.SpectrumAnalyzer();
SA.NumInputPorts = 1;
SA.SampleRate = fs;
SA.SpectrumType = 'Power';
SA.PowerUnits = 'dBW';
SA.Method = 'Filter Bank';
SA.PlotAsTwoSidedSpectrum = true;
SA.FrequencyScale = 'Linear';
% SA.RBWSource = 'Property';
% SA.RBW = 100;
% SA.FrequencySpan = 'Start and stop frequencies';
% SA.StartFrequency = 0.5*fsig;
% SA.StopFrequency = fs/2;
SA.YLimits = [-200,50];
SA(x4.')