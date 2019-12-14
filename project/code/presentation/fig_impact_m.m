% fig_compareFilters.m

close all
clear

M = 1024;           % Number of FFT points
order = 2*M*2-1;
% TMUX parameter definition


% Prototype filter definition
% load(sprintf('prototype/prototype (M = %d, m = %d).mat',numFFT,m))
load(sprintf('prototype/prototype from Kaiser (M = %d, m = 1).mat',M))
filt1 = p0(200,:);
load(sprintf('prototype/prototype from Kaiser (M = %d, m = 2).mat',M))
filt2 = p0(200,:);


flist = ["Rect","Kaiser m=1, order=2047", "Kaiser m=2, order=4095" ];

p0 = p0/sqrt(2*M); % The prototype filter

%%--plot spectrum for prototype filter

hfvt = fvtool(rectwin(M)./sum(rectwin(M)));

addfilter(hfvt, filt1./sum(filt1))
addfilter(hfvt, filt2./sum(filt2))

% titleStr = sprintf('Filter Comparison (M = %d, m = %d)',M,m);
% title(titleStr);
xticks([0:1/(M/2):1])
for i = 1:M/2
    k = i-1;
    names(i) = sprintf("%d*2\\pi/M",k);
end
names = cellstr(names);
set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend( flist);
set(gca, 'fontsize',13)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
xlim([0,2*2*pi/M])
% set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 6 4])



%% Export the prototype filter
savefile=0;
if(savefile)
saveas(gcf, sprintf('figures/Effectm.pdf'));
!pdfcrop figures/Effectm.pdf figures/Effectm.pdf
end

