% Build the prototype filter based on the Kaiser window
clc
clf
close all
clearvars

%% Prototype filter parameter definition
% M = 16; % The number of channels (for demonstration)
M = 1024; % The number of channels
m = 2; % The length of each polyphase component
N = 2*m*M-1; % The order of the prototype filter
L = N+1; % The length of the prototype filter

%% Brute-force search
% Choose alpha (sidelobe attenuation) such that the power complementary condition is satisfied to the maximum extent
% alpha = [1e3,7e3];
% for i = 1:length(alpha)
%     diff(i) = diffFromPowComp2(alpha(i),M);
% end
% [minResult,minidx] = min(diff);
% alphaOpt = alpha(minidx);
% alphaOpt = 85; % (for demonstration)
alphaOpt = 250;

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
% set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend(hfvt, 'The prototype filter', sprintf('The Kaiser window (\\alpha = %d)',alphaOpt), 'The rectangular window');

%% Export the prototype filter
save(sprintf("prototype (M = %d, m = %d).mat",M,m),'p0');
saveas(gca,[titleStr,'.jpg'])