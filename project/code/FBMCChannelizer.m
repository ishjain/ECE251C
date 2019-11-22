% An FBMC channelizer in direct, polyphase, and overlap-add form
clc
clf
close all
clearvars

M = 4; % The number of channels
% L = 2; % The number of frames
xSerr = [1, 2, 3, 4, 5, 6, 7, 8, 9];
xSeri = [31, 55, 1.2, 4.5, 66, 78, 45, 13, 11];
xSer = complex(xSerr,xSeri); % Received samples from the channel
W = exp(-1j*2*pi/M);
h0 = [1,2,3,4,3,5,1]; % Some arbitrary prototype filter

%% A channelizer in direct form
% Loop over each channel
for idx = 1:M
    k = idx-1;
    h(idx,:) = h0.*W.^(-k*[0:length(h0)-1]); % The filter
    r(idx,:) = conv(h(idx,:),xSer); % The filtered signal at each channel
    sHat(idx,:) = downsample(r(idx,:),M); % The received symbols
end

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 2 polyphase is used)
% The polyphase components of the prototype filter go to each channel
numRow = M;
numCol = ceil(length(h0)/M);
R = reshape([h0,zeros(1,numRow*numCol-length(h0))],numRow,numCol); 
R = flipud(R); % The type 2 polyphase component of the prototype filter
FF = dftmtx(M); % The DFT matrix

% S/P conversion
L = ceil((length(xSer)-1)/M);
xPar2 = reshape([xSer(2:end),zeros(1,M*L-(length(xSer)-1))],M,L);
xPar2 = [[zeros(M-1,1);xSer(1)],xPar2];

% Loop over each channel for filtering
for idx = 1:M
    sHat2(idx,:) = conv(R(idx,:),xPar2(idx,:));
end

% Shuffling
SM = [zeros(1,M-1),1;eye(M-1),zeros(M-1,1)]; % The shuffling matrix
sHat2 = SM*sHat2;

% Compute the DFT
sHat2 = FF*sHat2;

%% A channel synthesizer in polyphase form (The prototype filter is upconverted and type 1 polyphase is used)
% The polyphase components of the prototype filter go to each channel
numRow = M;
numCol = ceil(length(h0)/M);
E = reshape([h0,zeros(1,numRow*numCol-length(h0))],numRow,numCol); 
FF = dftmtx(M);
IFF = conj(FF); % The IDFT matrix

% S/P conversion
L = ceil((length(xSer)-1)/M);
xPar1 = reshape([xSer(2:end),zeros(1,M*L-(length(xSer)-1))],M,L);
xPar1 = flipud(xPar1);
xPar1 = [[xSer(1);zeros(M-1,1)],xPar1];

% Loop over each channel for filtering
for idx = 1:M
    sHat1(idx,:) = conv(xPar1(idx,:),E(idx,:));
end

% Compute the IDFT
sHat1 = IFF*sHat1;

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