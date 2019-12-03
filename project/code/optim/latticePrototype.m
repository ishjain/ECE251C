% From angles to prototype
clc
clf
close all
clearvars

%% Prototype filter parameter definition
M = 16; % M channels for the FB, and 2M polyphase components for the prototype filter
m = 1; % m taps for each polyphase component
N = 2*m*M-1; % The order of the prototype filter

% theta = pi/4*ones(M/2,m); % The parameters to optimize, one row per channel
theta = pi/4*ones(M/2,1);
G(1:M/2) = cos(theta);
if isrow(G)
    G = G.';
end
G(M/2+1:M) = flipud(sin(theta));
G(M+1:M+M/2) = fliplr(flipud(G(M/2+1:M)));
G(M+M/2+1:2*M) = fliplr(flipud(G(1:M/2)));
p = G; % The prototype filter

% Compute the DFT
nfft = 1024*M;
x_ax = 0:1/nfft:1-1/nfft;
P = mag2db(abs(fft(p,nfft)));
plot(x_ax,P)