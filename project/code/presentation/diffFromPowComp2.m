function val = diffFromPowComp2(alpha,M)
% m = 2
L = 2*2*M; % The length of the prototype filter
if (alpha >= 21) && (alpha <= 50)
    beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
elseif alpha > 50
    beta = 0.1102*(alpha-8.7);
end
kw = kaiser(L,beta); % A Kaiser window    
% Extract the kth and the (M+k)th polyphase components out of the 2M
% polyphase components
for i = 1:M/2
    k = i-1;
    gk(i,:) = downsample(kw,2*M,k);
    gMplusk(i,:) = downsample(kw,2*M,M+k);
end
% The "cosines" and "sines" of the sums and differences of the angles
cd = gk(:,1)+gk(:,2);
sd = gMplusk(:,2)+gMplusk(:,1);
cs = gk(:,1)-gk(:,2);
ss = gMplusk(:,2)-gMplusk(:,1);
% The sum of squares of cosine and sine terms (ideally should be the same for each polyphase component)
% The cosines and sines of the sums and differences of the angles should
% also be ideally the same
sumSqSum = cs.^2+ss.^2;
sumSqDiff = cd.^2+sd.^2;
% Compute the average of the sum of squares for each polyphase component
avgSum = mean(sumSqSum);
avgDiff = mean(sumSqDiff);
% Compute the dB difference between these two terms
dBDiff = mag2db(avgSum) - mag2db(avgDiff);
val = abs(dBDiff);