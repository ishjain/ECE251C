clc
clf
close all
clearvars

M = 16;
m = 1;
diff = @(alpha) diffFromPowComp(alpha,M,m);
diff(40)