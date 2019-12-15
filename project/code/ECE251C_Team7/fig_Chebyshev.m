%% Chebyshev
% fig_compareFilters.m

close all
clear
M = 1024;           % Number of FFT points
m = 2;
N = 2*m*M-1;
load(sprintf('prototype/prototype from Chebyshev (M = %d, m = %d).mat',M,m))
p0 = p0/sqrt(2*M); % The prototype filter
%%--plot spectrum for prototype filter
idx1 = 1;
idx2 = 40;
idx3 = 81;
hfvt = fvtool(p0(idx1,:)./sum(p0(idx1,:)));
addfilter(hfvt, p0(idx2,:)./sum(p0(idx2,:)))
addfilter(hfvt, p0(idx3,:)./sum(p0(idx3,:)))
titleStr = sprintf('Chebyshev Prototype Filter (M = %d, m = %d)',M,m);
title(titleStr);
xticks([0:1/(M/2):1])
for i = 1:M/2
    k = i-1;
    names(i) = sprintf("%d*2\\pi/M",k);
end
names = cellstr(names);
set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend(hfvt, sprintf('The prototype filter r = %d',r(idx1)), sprintf('The prototype filter r = %d',r(idx2)), sprintf('The prototype filter r = %d',r(idx3)));
set(gca, 'fontsize',13)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
xlim([0,2*2*pi/M])
%% Export the prototype filter
savefile=0;
if(savefile)
saveas(gcf, sprintf('figures/ChebyshevPrototypeFilter_Effectr.pdf',M,m));
!pdfcrop figures/ChebyshevPrototypeFilter_Effectr.pdf figures/ChebyshevPrototypeFilter_Effectr.pdf
end