% Build the prototype filter based on the Kaiser window
clc
clf
close all
clearvars

%% Prototype filter parameter definition
M = 1024; % The number of channels
m = 1; % The length of each polyphase component
N = 2*m*M-1; % The order of the prototype filter
L = N+1; % The length of the prototype filter

%% Brute-force search
% Choose alpha (sidelobe attenuation) such that the power complementary condition is satisfied to the maximum extent
alpha = 22:0.01:80;
for i = 1:length(alpha)
    diff(i) = diffFromPowComp1(alpha(i),M);
end
[minResult,minidx] = min(diff);
alphaOpt = alpha(minidx);

%% Compare the prototype filter with the Kaiser window and rectangular window
[kw,p0,theta] = estAngle(alphaOpt,M,m);
hfvt = fvtool(p0./sum(p0));
addfilter(hfvt, kw./sum(kw));
addfilter(hfvt, rectwin(M)./sum(rectwin(M)))
titleStr = sprintf('Filter Comparison (M = %d, m = %d)',M,m);
title(titleStr);
% xticks([0:1/(M/2):1])
for i = 1:M/2
    k = i-1;
    names(i) = sprintf("%d*2\\pi/M",k);
end
names = cellstr(names);
set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend(hfvt, 'The prototype filter', 'The Kaiser window', 'The rectangular window');
xlim([0,2*2*pi/M])
%% Export the prototype filter
savefile=0;
if(savefile)
save(sprintf("prototype (M = %d, m = %d).mat",M,m),'p0');
saveas(gca,[titleStr,'.jpg'])
end