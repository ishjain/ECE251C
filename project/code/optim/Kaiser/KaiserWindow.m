clc
clf
close all
clearvars

M = 16;
m = 1;
N = 2*m*M-1;
L = N+1;
alpha = 35.76;
beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
p0 = kaiser(L,beta);
g = conv(p0,p0);

% fvtool(p0)
% fvtool(g)

c = p0(1:M/2);
s = flipud(p0(M/2+1:M));

c.^2+s.^2

p0 = p0/sum(p0);
fvtool(p0)

% Return the approximate angle values
theta = atan2(s,c);

% Build a prototype filter out of the angles
G(1:M/2) = cos(theta);
if isrow(G)
    G = G.';
end
G(M/2+1:M) = flipud(sin(theta));
G(M+1:M+M/2) = fliplr(flipud(G(M/2+1:M)));
G(M+M/2+1:2*M) = fliplr(flipud(G(1:M/2)));
p = G; % The prototype filter

p = p/sum(p);
fvtool(p)