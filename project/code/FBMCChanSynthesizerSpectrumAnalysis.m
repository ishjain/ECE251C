% An FBMC channel synthesizer in various forms and the spectrum of the
% output
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
Kov = 4;
n = 0:Kov*M-1;
b1 = 0.97195983;
b2 = 1/sqrt(2);
b3 = 0.23514695;
f0 = 1+2*(b1*cos(2*pi*n/(Kov*M))+...
    b2*cos(2*2*pi*n/(Kov*M))+...
    b3*cos(3*2*pi*n/(Kov*M)));

%% A channel synthesizer in direct form
% Loop over each channel
for idx = 1:M
    k = idx-1;
    f(idx,:) = f0.*W.^(-k*[0:length(f0)-1]); % The filter
    w(idx,:) = upsample(S(idx,:),M); % The symbols after upsampling
    p(idx,:) = conv(w(idx,:),f(idx,:)); % The filtered symbols
end
x = sum(p); % Sum over all channels

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 2 polyphase is used)
% % The polyphase components of the prototype filter go to each channel
% numRow = M;
% numCol = ceil(length(f0)/M);
% R = reshape([f0,zeros(1,numRow*numCol-length(f0))],numRow,numCol); 
% R = flipud(R); % The type 2 polyphase component of the prototype filter
% FF = dftmtx(M); % The DFT matrix
% 
% SS = FF*S;
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
% % P/S conversion
% x2 = flipud(Q);
% x2 = x2(:).';

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 1 polyphase is used)
% % This is the approach appearing in comm textbooks
% % The polyphase components of the prototype filter go to each channel
% numRow = M;
% numCol = ceil(length(f0)/M);
% R = reshape([f0,zeros(1,numRow*numCol-length(f0))],numRow,numCol); 
% FF = dftmtx(M);
% IFF = conj(FF); % The IDFT matrix
% 
% SS = IFF*S;
% 
% % Loop over each channel for filtering
% for idx = 1:M
%     k = idx-1;
%     Q(idx,:) = conv(R(idx,:),SS(idx,:));
% end
% 
% % P/S conversion
% x3 = Q(:).';

%% Overlap-add in the time domain
% FF = conj(dftmtx(M)); % The IDFT matrix
% repnum = ceil(length(f0)/M);
% Q = FF*S;
% Q = repmat(Q,repnum,1);
% Q = Q(1:length(f0),:);
% Q = Q.*f0.';
% for idx = 1:size(Q,2)
%     i = idx-1;
%     x4(idx,1:i*M+length(f0)) = [zeros(1,i*M),Q(:,idx).'];
% end
% x4 = sum(x4);

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
SA(x.')