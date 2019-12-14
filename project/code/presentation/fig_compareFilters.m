% fig_compareFilters.m

close all

M = 1024;           % Number of FFT points

% TMUX parameter definition
m = 2;
N = 2*m*M-1;

% Prototype filter definition
% load(sprintf('prototype/prototype (M = %d, m = %d).mat',numFFT,m))
load(sprintf('prototype/prototype from Kaiser (M = %d, m = %d).mat',M,m))
filt(1,:) = p0(200,:);
load(sprintf('prototype/prototype from Chebyshev (M = %d, m = %d).mat',M,m))
filt(2,:) = p0(40,:);

load(sprintf('prototype/Misc Prototype 1 (M = %d, m = %d).mat',M,m))
flist = ["Rect","Kaiser", "Chebyshev", "barthannwin","blackman","blackmanharris","bohmanwin","gausswin","flattopwin","hamming","hann","nuttallwin","parzenwin"];

goodlist=[1,2,3,10, 5, 6];

for ii=1:10
filt(2+ii,:) = p0(ii,:);
end

% p0=p0(200,:);
p0 = p0/sqrt(2*M); % The prototype filter

%%--plot spectrum for prototype filter

hfvt = fvtool(rectwin(M)./sum(rectwin(M)));

for  ind  =1:length(goodlist)
    ii=goodlist(ind);
addfilter(hfvt, filt(ii,:)./sum(filt(ii,:)))
end

titleStr = sprintf('Filter Comparison (M = %d, m = %d)',M,m);
title(titleStr);
% xticks([0:1/(M/2):1])
% for i = 1:M/2
%     k = i-1;
%     names(i) = sprintf("%d*2\\pi/M",k);
% end
% names = cellstr(names);
% set(gca,'xtick',[0:1/(M/2):1],'xticklabel',names)
xlabel('Normalized Frequency (rad/sample)')
legend( flist(goodlist), 'location', 'bestoutside');
set(gca, 'fontsize',13)
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
xlim([0,2*2*pi/M])
set(gcf,'PaperUnits', 'inches', 'paperposition', [0 0 10 5])


% set(gcf,'Units','Inches');
% 
 print(gcf,'filename','-dpdf','-r0')

%% Export the prototype filter
savefile=0;
if(savefile)
saveas(gcf, sprintf('figures/Comparefilters.pdf'));
!pdfcrop figures/Comparefilters.pdf figures/Comparefilters.pdf
end


