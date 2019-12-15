% fig_compareFilters.m

close all
clear
M = 1024;           % Number of FFT points
m = 2;
N = 2*m*M-1;

load(sprintf('prototype/prototype from Kaiser (M = %d, m = %d).mat',M,m))

p0 = p0/sqrt(2*M); % The prototype filter

%%--plot spectrum for prototype filter
idx1 = 1;
idx2 = 200;
idx3 = 271;
hfvt = fvtool(p0(idx1,:)./sum(p0(idx1,:)));
addfilter(hfvt, p0(idx2,:)./sum(p0(idx2,:)))
addfilter(hfvt, p0(idx3,:)./sum(p0(idx3,:)))
titleStr = sprintf('Kaiser Prototype Filter (M = %d, m = %d)',M,m);
title(titleStr);
xticks([0:1/(M/2):1])
for i = 1:M/2
    k = i-1;
    names(i) = sprintf("%d*2\\pi/M",k);
end
names = cellstr(names);
set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend(hfvt, sprintf('The prototype filter \\alpha = %d',alpha(idx1)), sprintf('The prototype filter \\alpha = %d',alpha(idx2)), sprintf('The prototype filter \\alpha = %d',alpha(idx3)));
set(gca, 'fontsize',13)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
xlim([0,2*2*pi/M])
%% Export the prototype filter
savefile=0;
if(savefile)
saveas(gcf, sprintf('figures/KaiserPrototypeFilter_EffectAlpha.pdf',M,m));
!pdfcrop figures/KaiserPrototypeFilter_EffectAlpha.pdf figures/KaiserPrototypeFilter_EffectAlpha.pdf
end


