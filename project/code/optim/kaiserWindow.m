clc
clf
close all
clearvars

M = 16;
m = 1;
N = 2*m*M-1;
L = N+1;
alpha = 40;
beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
p0 = kaiser(L,beta);
g = conv(p0,p0);

% fvtool(p0)
% fvtool(g)

c = p0(1:M/2);
s = flipud(p0(M/2+1:M));

c.^2+s.^2