function val = PStopMax(theta)
M = 16;
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
range = nfft/(4*M):nfft/2;
P = mag2db(abs(fft(p,nfft)));
maxInRange = max(P(range));
val = maxInRange-P(1);