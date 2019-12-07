function val = diffFromPowComp1(alpha,M)
% m = 1
L = 2*M; % The length of the prototype filter
if (alpha >= 21) && (alpha <= 50)
    beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
elseif alpha > 50
    beta = 0.1102*(alpha-8.7);
end
kw = kaiser(L,beta); % A Kaiser window
c = kw(1:M/2); % The "cosine" terms
s = flipud(kw(M/2+1:M)); % The "sine" terms
sumSq = c.^2+s.^2; % The sum of squares of cosine and sine terms (ideally should be the same for each polyphase component)
sumSqdB = mag2db(sumSq);
dBDiff = max(sumSqdB)-sumSqdB; % Find the dB difference between the maximum of sum of squares and each of the sum of squares
val = max(dBDiff); % The maximal value of the above difference