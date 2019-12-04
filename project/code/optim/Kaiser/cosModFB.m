% An FB from the Kaiser window prototype
clc
clf
close all
clearvars

%% Filter bank parameter definition
M = 16;
m = 1;
N = 2*m*M-1;

%% Prototype filter definition
thetaOptim = [1.34835848234932;1.28727140093828;1.22112655920844;1.15007561376553;1.07448741731735;0.994982113744842;0.912438788862042;0.827965788036032];
theta = thetaOptim;
G(1:M/2) = cos(theta);
if isrow(G)
    G = G.';
end
G(M/2+1:M) = flipud(sin(theta));
G(M+1:M+M/2) = fliplr(flipud(G(M/2+1:M)));
G(M+M/2+1:2*M) = fliplr(flipud(G(1:M/2)));
p0 = G;
p0 = p0/sum(p0); % The prototype filter
if iscolumn(p0)
    p0 = p0.';
end

%% Filter bank definition
n = 0:N;
for idx = 1:M
    k = idx-1;
    thetak = (-1)^k*pi/4;
    h(idx,:) = 2*p0.*cos(pi/M*(k+1/2)*(n-N/2)+thetak);
    f(idx,:) = 2*p0.*cos(pi/M*(k+1/2)*(n-N/2)-thetak);
end

%% Test for PR
x = rand(1,100);
for idx = 1:M
    hout(idx,:) = conv(x,h(idx,:));
    DSout(idx,:) = downsample(hout(idx,:),M);
    USout(idx,:) = upsample(DSout(idx,:),M);
    fout(idx,:) = conv(USout(idx,:),f(idx,:));
end
xHat = sum(fout);
xHat = xHat*x(1)/xHat(32);