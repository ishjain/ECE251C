% Check the orthogonality condition of prototype filters
clc
clf
close all
clearvars

M = 1024; % The number of channels
F = 15e3; % Subcarrier Spacing: 15kHz
fs = M*F; % Sampling Rate: 1024*15kHz

dt = 1/fs; % Sampling interval

%% Define Prototype Filters! Orthogonal for T=T0 and F=2/T0!
% PHYDYAS pulse for overlapping factor 4
p_PHYDYAS_O4 = @(t,T0) ((t<=(4*T0/2))&(t>-4*T0/2)).*(1+2*(...
    0.97195983 * cos(2*pi*1/4*t/T0) + ...
    sqrt(2)/2  * cos(2*pi*2/4*t/T0) + ...
    0.23514695 * cos(2*pi*3/4*t/T0) ...
    ))/sqrt(4^2*T0);

% PHYDYAS pulse for overlapping factor 8
% p_PHYDYAS_O8 = @(t,T0) ((t<=(8*T0/2))&(t>-8*T0/2)).*(1+2*(...
%     0.99932588 * cos(2*pi*1/8*t/T0) + ...
%     0.98203168 * cos(2*pi*2/8*t/T0) + ...
%     0.89425129 * cos(2*pi*3/8*t/T0) + ...
%     sqrt(2)/2  * cos(2*pi*4/8*t/T0) + ...
%     0.44756522 * cos(2*pi*5/8*t/T0) + ...
%     0.18871614 7.934249631012864e-08 - 8.208809704409487e-18i* cos(2*pi*6/8*t/T0) + ...
%     0.03671221 * cos(2*pi*7/8*t/T0) ...
%     ))/sqrt(8^2*T0);

% Root raised cosine filter in time (low latency)
p_timeRRC = @(t,T0) ((t<=(T0/2))&(t>-T0/2)).*sqrt(1+(...
      cos(pi*1*2*t/T0)  ...
     ))/sqrt(T0);
 
% Hermite prototype filter for overlapping factor 8
O  = 8; 
p_Hermite = @(t,T0) ((t<=(O*T0/2))&(t>-O*T0/2)).* ...
    1/sqrt(T0).*exp(-pi*(t./(T0/sqrt(2))).^2) .* (...
    1.412692577 + ...
    -3.0145e-3 .*...
                ((12+(-48).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+16.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4 ) )+ ...
    -8.8041e-6 .*...
                (1680+(-13440).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+13440.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-3584).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+256.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8 )+ ...
    -2.2611e-9  .*... 
                 (665280+(-7983360).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+13305600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-7096320).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+1520640.* ...
                    (sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-135168).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+4096.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12 )+ ...
    -4.4570e-15 .*... 
                 (518918400+(-8302694400).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+19372953600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-15498362880).* ...
                   (sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+5535129600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-984023040).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+89456640.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12+( ...
                   -3932160).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^14+65536.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^16 )+ ...
     1.8633e-16 .*...
                 (670442572800+(-13408851456000).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+40226554368000.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+( ...
                   -42908324659200).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+21454162329600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-5721109954560).* ...
                   (sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+866834841600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12+(-76205260800).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^14+3810263040.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^16+ ...
                   (-99614720).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^18+1048576.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^20 ));
                
% Reference rectangular prototype filter, orthogonal for T=T0, F=1/T0!
p_Rectangular = @(t,T0) (1/sqrt(T0).*((t<=(T0/2))&(t>-T0/2)));

%% Ambiguity function definition
p_PulseDelayMod = @(t,T0,l,k,F) p_PHYDYAS_O4(t-k*T0,T0).*exp(1j*2*pi*l*F*(t-k*T0));
p_Ambiguity = @(t,T0,k1,k2,l1,l2,F) p_PulseDelayMod(t,T0,l1,k1,F).*conj(p_PulseDelayMod(t,T0,l2,k2,F));

%% Plot Filters
% T0  = 1/F; % Frame interval
% t   = (-200*T0/2:dt:(200*T0/2)).'; t(end)=[]; % so that fft becomes real
% 
% % Time domain
% % p_PHYDYAS_O4_Samples    = p_PHYDYAS_O4(t,T0);
% % p_PHYDYAS_O8_Samples    = p_PHYDYAS_O8(t,T0);
% % p_timeRRC_Samples       = p_timeRRC(t,T0);
% % p_Hermite_Samples       = p_Hermite(t,T0);
% p_Rectangular_Samples   = p_Rectangular(t,T0);
% % p_RectCaller_Samples   = p_RectCaller(t,T0);
% p_Ambiguity_Samples = p_Ambiguity(t,T0,1,0,2,0,F);
% 
% % Plot the prototype filters in time
% figure();
% % plot(t/T0,p_PHYDYAS_O8_Samples*sqrt(T0),'red');
% % hold on;
% % plot(t/T0,p_timeRRC_Samples*sqrt(T0),'magenta');
% % plot(t/T0,p_Hermite_Samples*sqrt(T0),'blue');
% plot(t/T0,p_Rectangular_Samples*sqrt(T0),'black');
% hold on;
% % plot(t/T0,p_RectCaller_Samples*sqrt(T0),'red');
% plot(t/T0,p_Ambiguity_Samples*sqrt(T0),'yellow');
% xlim([-3 3]);
% ylabel('p(t)');
% xlabel('Normalized Time, t/T0');
% legend({'PHYDYAS','RRC','Hermite','Rectangular'});

%% Check orthogonality
T0  = 1/F; % Frame interval
l1 = 0;
l2 = 0;
k1 = 0;
l2list=0:10;
k2list = 0:10;
% k2 = 0;
for k2i = 1:length(k2list)
    for l2i=1:length(l2list)
        l2 = l2list(l2i);
        k2=k2list(k2i);
        p_OrFcn = @(t) p_Ambiguity(t,T0,k1,k2,l1,l2,F);
        loBound = min(k1,k2)*T0-T0*2;
        upBound = max(k1,k2)*T0+T0*2;
        Result(k2i,l2i) = integral(p_OrFcn,loBound,upBound);
    end
    % fprintf('Result = %.8f\n',Result)
end
figure;
[X,Y]= meshgrid(k2list,l2list);
figure;
mesh(X,Y,abs(Result));
figure;
imshow(abs(Result))
figure;
stem(abs(Result))
max(abs(Result))