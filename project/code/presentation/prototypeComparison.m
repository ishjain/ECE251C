%% Compare the two prototype filter with the Kaiser window and rectangular window
clc
clf
clearvars
close all

M = 16;
m1 = 1;
m2 = 2;
alpha1 = 35.76;
alpha2 = 120;
[kw1,p1,theta1] = estAngle(alpha1,M,m1);
[kw2,p2,theta2] = estAngle(alpha2,M,m2);
hfvt = fvtool(p1./sum(p1));
addfilter(hfvt, p2./sum(p2));
addfilter(hfvt, rectwin(M)./sum(rectwin(M)))
titleStr = sprintf('Filter Comparison (M = %d)',M);
title(titleStr);

for i = 1:M/2
    k = i-1;
    names(i) = sprintf("%d/M",k);
end
names = cellstr(names);
set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (f/f_s)')
legend(hfvt, sprintf('The prototype filter (m = %d, \\alpha = %.2f)',m1,alpha1), sprintf('The prototype filter (m = %d, \\alpha = %d)',m2,alpha2), 'The rectangular window');