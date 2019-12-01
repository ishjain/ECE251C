% The PHYDYAS prototype filter in discrete time
clc
clf
close all
clearvars

M = 1024;
h = @(n) ((n>=0)&(n<M)).*1/sqrt(M);
f = @(n) h(-n);
hprod = @(m,n) h(m).*f(n-m);

%% Check if the product of the prototype is Nyquist or not
% n = -M:M;
% mm = -M:M;
% for i = 1:length(n)
%     gg = @(m) hprod(m,n(i));
%     g(i) = sum(gg(mm));
% end
% 
% g = g/max(g);
% idx1 = find(n==-M);
% idx2 = find(n==M);
% idx3 = find(n==0);
% 
% figure()
% plot(n,g)
% hold on
% plot(n(idx1),g(idx1),'ro')
% plot(n(idx2),g(idx2),'ro')
% plot(n(idx3),g(idx3),'ro')
% text(n(idx1),g(idx1)+0.1,sprintf('g[%d] = %d',n(idx1),g(idx1)))
% text(n(idx2),g(idx2)+0.1,sprintf('g[%d] = %d',n(idx2),g(idx2)))
% text(n(idx3),g(idx3),sprintf('g[%d] = %d',n(idx3),g(idx3)))
% hold off
% xlabel('Sample index')
% ylabel('Amplitude')

%% Now if the filters are derived by exponential modulation
k1 = 1;
k2 = 100;
hk1 = @(n,k1) h(n).*exp(1j*k1*2*pi/M*n);
fk2 = @(n,k2) f(n).*exp(1j*k2*2*pi/M*n);
hk1k2prod = @(n,m,k1,k2) hk1(m,k1).*fk2(n-m,k2);

n = -M:M;
mm = -M:M;
for i = 1:length(n)
    ggk1k2 = @(m) hk1k2prod(n(i),m,k1,k2);
    gk1k2(i) = sum(ggk1k2(mm));
end

gk1k2 = gk1k2/max(gk1k2);
idx1 = find(n==-M);
idx2 = find(n==M);
idx3 = find(n==0);

figure()
plot(n,abs(gk1k2))
hold on
plot(n(idx1),abs(gk1k2(idx1)),'ro')
plot(n(idx2),abs(gk1k2(idx2)),'ro')
plot(n(idx3),abs(gk1k2(idx3)),'ro')
text(n(idx1),abs(gk1k2(idx1))+0.1,sprintf('g[%d] = %d',n(idx1),abs(gk1k2(idx1))))
text(n(idx2),abs(gk1k2(idx2))+0.1,sprintf('g[%d] = %d',n(idx2),abs(gk1k2(idx2))))
text(n(idx3),abs(gk1k2(idx3))+0.1,sprintf('g[%d] = %d',n(idx3),abs(gk1k2(idx3))))
hold off
xlabel('Sample index')
ylabel('Amplitude')

%% What about cosine modulation
% N = K*M-1;
% k1 = 5;
% k2 = 5;
% hk1 = @(n,k1) 2*h(n).*cos(pi/M*(k1+1/2)*(n-N/2)+(-1)^k1*pi/4);
% fk2 = @(n,k2) 2*f(n).*cos(pi/M*(k2+1/2)*(n-N/2)-(-1)^k2*pi/4);
% hk1k2prod = @(n,m,k1,k2) hk1(m,k1).*fk2(n-m,k2);
% 
% n = -K*M:K*M;
% mm = -K*M:K*M;
% for i = 1:length(n)
%     ggk1k2 = @(m) hk1k2prod(n(i),m,k1,k2);
%     gk1k2(i) = sum(ggk1k2(mm));
% end
% 
% gk1k2 = gk1k2/max(gk1k2);
% idx1 = find(n==-M);
% idx2 = find(n==M);
% idx3 = find(n==-2*M);
% idx4 = find(n==2*M);
% idx5 = find(n==0);
% 
% figure()
% plot(n,abs(gk1k2))
% hold on
% plot(n(idx1),abs(gk1k2(idx1)),'ro')
% plot(n(idx2),abs(gk1k2(idx2)),'ro')
% plot(n(idx3),abs(gk1k2(idx3)),'ro')
% plot(n(idx4),abs(gk1k2(idx4)),'ro')
% plot(n(idx5),abs(gk1k2(idx5)),'ro')
% text(n(idx1),gk1k2(idx1)+0.1,sprintf('g_{k1k2}[%d] = %d',n(idx1),gk1k2(idx1)))
% text(n(idx2),gk1k2(idx2),sprintf('gk1k2[%d] = %d',n(idx2),gk1k2(idx2)))
% text(n(idx3),gk1k2(idx3)+0.2,sprintf('g1k2[%d] = %d',n(idx3),gk1k2(idx3)))
% text(n(idx4),gk1k2(idx4)+0.1,sprintf('g1k2[%d] = %d',n(idx4),gk1k2(idx4)))
% text(n(idx5),gk1k2(idx5),sprintf('g1k2[%d] = %d',n(idx5),gk1k2(idx5)))
% hold off
% xlabel('Sample index')
% ylabel('Amplitude')
% 
% figure()
% plot(mag2db(abs(fft(hk1(n,k1)))))