% m = 2
clc
clf
close all
clearvars

%% Window to angle
M = 16;
m = 2;
N = 2*m*M-1;
L = N+1;
% 
% idx = 1;
% for alpha = 100:1:140
%     if alpha > 50
%         beta = 0.1102*(alpha-8.7);
%     elseif alpha >= 21
%         beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
%     end
%     
%     p0 = kaiser(L,beta);
%     p0 = p0/sum(p0);
%     % fvtool(p0)
% 
%     for i = 1:M/2
%         k = i-1;
%         gk(i,:) = downsample(p0,2*M,k);
%         gMplusk(i,:) = downsample(p0,2*M,M+k);
%     end
% 
%     cm(:,idx) = sum(gk,2);
%     sm(:,idx) = sum(gMplusk,2);
%     cp(:,idx) = gk(:,1)-gk(:,2);
%     sp(:,idx) = gMplusk(:,2)-gMplusk(:,1);
%     psqsum(:,idx) = cp(:,idx).^2+sp(:,idx).^2;
%     msqsum(:,idx) = cm(:,idx).^2+sm(:,idx).^2;
%     idx = idx+1;
% end
% pmean = mean(psqsum);
% mmean = mean(msqsum);
% diff = pmean-mmean;

%% Try this particular alpha
alpha = 100;
beta = 0.1102*(alpha-8.7);   
p0 = kaiser(L,beta);
p0 = p0/sum(p0);
fvtool(p0)

for i = 1:M/2
    k = i-1;
    gk(i,:) = downsample(p0,2*M,k);
    gMplusk(i,:) = downsample(p0,2*M,M+k);
end

cm = sum(gk,2);
sm = sum(gMplusk,2);
cp = gk(:,1)-gk(:,2);
sp = gMplusk(:,2)-gMplusk(:,1);
psqsum = cp.^2+sp.^2;
msqsum = cm.^2+sm.^2;


% msqsumdB = pow2db(msqsum);
% diff = max(max(sqsumdB)-msqsumdB);
% gk = gk.*sqrt(1./msqsum);
% gMplusk = gMplusk.*sqrt(1./msqsum);
% cm = sum(gk,2);
% sm = sum(gMplusk,2);
% sqsum = cm.^2+sm.^2;
% cm = cm.*sqrt(1./sqsum);
% sm = sm.*sqrt(1./sqsum);
% msqsumdB = pow2db(msqsum);
% diff = max(max(sqsumdB)-msqsumdB);

thetam = atan2(sm,cm);
thetap = atan2(sp,cp);
theta0 = ((thetap+thetam)/2);
theta1 = ((thetap-thetam)/2);

%% Angle to polyphase
for idx = 1:M/2
    gg(idx,:) = [cos(theta0(idx))*cos(theta1(idx)), sin(theta0(idx))*sin(theta1(idx))];
    ggMplusk(idx,:) = [-cos(theta0(idx))*sin(theta1(idx)), sin(theta0(idx))*cos(theta1(idx))];
    gg(2*M+1-idx,:) = fliplr(gg(idx,:));
    gg(M+idx,:) = [-cos(theta0(idx))*sin(theta1(idx)), sin(theta0(idx))*cos(theta1(idx))];
    gg(M+1-idx,:) = fliplr(gg(M+idx,:));
end

% Polyphase to filter
pp0 = gg(:);
pp0 = pp0/sum(pp0);

fvtool(pp0)