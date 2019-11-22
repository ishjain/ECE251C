% An FBMC channel synthesizer in direct, polyphase, and overlap-add form
clc
clf
close all
clearvars

M = 4; % The number of channels
L = 2; % The number of frames
Sr = [1,2;...
    3,4;...
    5,6;...
    7,8];
Si = [31,55;...
    1.2,4.5;...
    66,78;...
    45,13];
S = complex(Sr,Si); % Symbols to transmit, row for channel and column for frame
W = exp(-1j*2*pi/M);
f0 = [1,2,3,4,3,5,1]; % Some arbitrary prototype filter

%% A channel synthesizer in direct form
% Loop over each channel
for idx = 1:M
    k = idx-1;
    f(idx,:) = f0.*W.^(-k*[0:length(f0)-1]); % The filter
    w(idx,:) = upsample(S(idx,:),M); % The symbols after upsampling
    p(idx,:) = conv(w(idx,:),f(idx,:)); % The filtered symbols
end
x1 = sum(p); % Sum over all channels

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 2 polyphase is used)
% The polyphase components of the prototype filter go to each channel
numRow = M;
numCol = ceil(length(f0)/M);
R = reshape([f0,zeros(1,numRow*numCol-length(f0))],numRow,numCol); 
R = flipud(R); % The type 2 polyphase component of the prototype filter
FF = dftmtx(M); % The DFT matrix

SS = FF*S;

% Shuffling
SM = [zeros(M-1,1),eye(M-1);1,zeros(1,M-1)];
SS = SM*SS;

% Loop over each channel for filtering
for idx = 1:M
    k = idx-1;
    Q(idx,:) = conv(R(idx,:),SS(idx,:));
end

% P/S conversion
x2 = flipud(Q);
x2 = x2(:).';

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 1 polyphase is used)
% This is the approach appearing in comm textbooks
% The polyphase components of the prototype filter go to each channel
numRow = M;
numCol = ceil(length(f0)/M);
R = reshape([f0,zeros(1,numRow*numCol-length(f0))],numRow,numCol); 
FF = dftmtx(M);
IFF = conj(FF); % The IDFT matrix

SS = IFF*S;

% Loop over each channel for filtering
for idx = 1:M
    k = idx-1;
    Q(idx,:) = conv(R(idx,:),SS(idx,:));
end

% P/S conversion
x3 = Q(:).';

%% Overlap-add in the time domain
FF = conj(dftmtx(M)); % The IDFT matrix
repnum = ceil(length(f0)/M);
Q = FF*S;
Q = repmat(Q,repnum,1);
Q = Q(1:length(f0),:);
Q = Q.*f0.';
for idx = 1:size(Q,2)
    i = idx-1;
    x4(idx,1:i*M+length(f0)) = [zeros(1,i*M),Q(:,idx).'];
end
x4 = sum(x4);