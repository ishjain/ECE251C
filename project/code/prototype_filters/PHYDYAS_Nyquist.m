% Check if PHYDYAS product filter satisfies Nyquist criterion
clc
clf
close all
clearvars

H1 = 0.97196;
H2 = 1/sqrt(2);
H3 = 0.235149;
K = 4;
T = 1;
h = @(t) ((t>=-K*T/2)&(t<K*T/2)).*...
    (1+2*(...
    H1*cos(1*2*pi/(K*T)*t)+...
    H2*cos(2*2*pi/(K*T)*t)+...
    H3*cos(3*2*pi/(K*T)*t)));
hprod = @(t,tau) h(tau).*h(t-tau);

%% Check if the product of the prototype is Nyquist or not
% t = -K*T:T/1e2:K*T;
% for i = 1:length(t)
%     gg = @(tau) hprod(t(i),tau);
%     g(i) = integral(gg,-K*T,K*T);
% end
% 
% idx1 = find(t==-T);
% idx2 = find(t==T);
% 
% figure()
% plot(t/T,g)
% hold on
% plot(t(idx1),g(idx1),'ro')
% plot(t(idx2),g(idx2),'ro')
% text(t(idx1),g(idx1)+1,sprintf('g(%d*T) = %d',t(idx1)/T,g(idx1)))
% text(t(idx2),g(idx2),sprintf('g(%d*T) = %d',t(idx2)/T,g(idx2)))
% hold off
% xlabel('Time (\timesT)')
% ylabel('Amplitude')

%% Now check across channels
F = 1/T;
k1 = 1;
k2 = 1;
hk = @(t,k) h(t).*exp(1j*k*2*pi*F*t);
hk1k2prod = @(t,tau,k1,k2) hk(tau,k1).*hk(t-tau,k2);

t = -K*T:T/1e2:K*T;
for i = 1:length(t)
    ggk1k2 = @(tau) hk1k2prod(t(i),tau,k1,k2);
    gk1k2(i) = integral(ggk1k2,-K*T,K*T);
end

figure()
plot(t/T,abs(gk1k2))
hold off
xlabel('Time (\timesT)')
ylabel('Amplitude')