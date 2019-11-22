% Comparing channel synthesizers in direct form and polyphase form
clc
clf
close all
clearvars

M = 4; % The number of channels
L = 2; % The number of frames
% S = complex(ones(M,L),zeros(M,L)); % Create a set of complex symbols to transmit
Sr = [1,2;...
    3,4;...
    5,6;...
    7,8];
Si = [31,55;...
    1.2,4.5;...
    66,78;...
    45,13];
S = complex(Sr,Si);
W = exp(-1j*2*pi/M);

%% A channel synthesizer in direct form
f0 = ones(1,M)/sqrt(M); % The prototype filter
% ff = conj(dftmtx(M))/sqrt(M);

% Loop over each channel
for idx = 1:M
    k = idx-1;
    f(idx,:) = f0.*W.^(-k*[0:M-1]); % The filter
    w(idx,:) = upsample(S(idx,:),M); % The symbols after upsampling
    p(idx,:) = conv(w(idx,:),f(idx,:)); % The filtered symbols
end
x1 = sum(p); % Sum over all channels

%% A channel synthesizer in polyphase form (version 1)
% The polyphase components of the prototype filter go to each channel
R = ones(M,1)/sqrt(M); 
FF2 = dftmtx(M); % The DFT matrix

% Compute DFT by looping over each frame
for idx = 1:L
    SS(:,idx) = FF2*S(:,idx); % Data symbols after DFT
end

% Shuffling
SM = [zeros(M-1,1),eye(M-1);1,zeros(1,M-1)];
SS = SM*SS;

% Loop over each channel for filtering
for idx = 1:M
    k = idx-1;
    Q(idx,:) = conv(R(idx,:),SS(idx,:));
end

% S/P conversion
x2 = flipud(Q);
x2 = x2(:).';

%% A channel synthesizer in polyphase form (version 2)
% The polyphase components of the prototype filter go to each channel
R = ones(M,1)/sqrt(M); 
FF = conj(dftmtx(M)); % The IDFT matrix

% Compute IDFT by looping over each frame
for idx = 1:L
    SS(:,idx) = FF*S(:,idx); % Data symbols after IDFT
end

% Loop over each channel for filtering
for idx = 1:M
    k = idx-1;
    Q(idx,:) = conv(R(idx,:),SS(idx,:));
end

% S/P conversion
x3 = flipud(Q);
x3 = x3(:).';