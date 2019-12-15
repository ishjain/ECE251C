% Compare different windows
clc
clf
close all
clearvars

%% Prototype filter parameter definition
M = 16; % The number of channels
% M = 1024; % The number of channels
m = 1; % The length of each polyphase component
N = 2*m*M-1; % The order of the prototype filter
L = N+1; % The length of the prototype filter

%% Generate prototype filters for each of the following windows
winType = ["barthannwin",... // 1
 "blackman",... // 2
 "blackmanharris",... // 3
 "bohmanwin",... // 4
 "gausswin",... // 5
 "flattopwin",... // 6
 "hamming",... // 7
 "hann",... // 8
 "nuttallwin",... // 9
 "parzenwin"];... // 10
i = 7;
w = eval(sprintf('%s(%d)',winType(i),L));
[p0,theta] = win2Prototype(w,M,m);
[f,D] = diffFromPowComp(w,M,m);
w = w*f;

%% Chebyshev window
r = 30;
wC = chebwin(L,r); % The Chebyshev window
[p0C,thetaC] = win2Prototype(wC,M,m);
[fC,DC] = diffFromPowComp(wC,M,m);
wC = wC*fC;

%% Kaiser window
alpha = 30;
if (alpha >= 21) && (alpha <= 50)
    beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
elseif alpha > 50
    beta = 0.1102*(alpha-8.7);
end
wK = kaiser(L,beta); % The Kaiser window
[p0K,thetaK] = win2Prototype(wK,M,m);
[fK,DK] = diffFromPowComp(wK,M,m);
wK = wK*fK;

%% Visualize the prototype filter
% windowType = "Hamming";
% windowType = "Chebyshev";
windowType = "Kaiser";
if strcmp(windowType,"Kaiser")
    pp = p0K;
    ww = wK;
elseif strcmp(windowType,"Chebyshev")
    pp = p0C;
    ww = wC;
else
    pp = p0;
    ww = w;
end
hfvt = fvtool(pp./sum(pp));
addfilter(hfvt, ww./sum(ww));
titleStr = sprintf('Magnitude Responses of %s Window Design Example \n (M = %d, m = %d)',windowType,M,m);
title(titleStr,'fontsize',10);
legend(hfvt, sprintf('The prototype filter from %s (p_0[n])',windowType),...
        sprintf('The original %s window (w[n])',windowType),'location','best');
xlim([0, 0.5])
set(gcf, 'WindowStyle','normal','Position',  [0, 0, 1000, 1000])
clc

%% Export the prototype filter
% save(sprintf("prototype (M = %d, m = %d).mat",M,m),'winType','p0','w');
% saveas(gca,[titleStr,'.jpg'])

%% txt
wktxt = [];
wMplusktxt = [];
thetaktxt = [];
Dktxt = [];
for i = 1:M/2
    wktxt = [wktxt,sprintf('%.4f & ',wK(i))];
    wMplusktxt = [wMplusktxt,sprintf('%.4f & ',wK(M+i))];
    thetaktxt = [thetaktxt,sprintf('%.4f & ',thetaK(i))];
    Dktxt = [Dktxt,sprintf('%.4f & ',DK(i))];
end
fprintf('%s\n',wktxt);
fprintf('%s\n',wMplusktxt);
fprintf('%s\n',thetaktxt);
fprintf('%s\n',Dktxt);