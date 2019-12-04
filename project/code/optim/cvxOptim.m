% From angles to prototype
clc
clf
close all
clearvars

%% Prototype filter parameter definition
M = 16; % M channels for the FB, and 2M polyphase components for the prototype filter
m = 1; % m taps for each polyphase component
N = 2*m*M-1; % The order of the prototype filter

%% Optimization definition
cvx_begin
variable theta(M/2)
minimize(cos(theta))
subject to
    min(theta)>=0
cvx_end


theta = optimvar('theta',M/2,1);
theta0.theta = pi/4*ones(M/2,1);
minVal = @(theta) min(theta);
maxVal = @(theta) max(theta);
objExpr = fcn2optimexpr(@PStopMax,theta);
minConstr = fcn2optimexpr(minVal,theta);
maxConstr = fcn2optimexpr(maxVal,theta);
convProb = optimproblem('Objective',objExpr,'Constraints',minConstr>=0);
[sol,fval,exitflag,output] = solve(convProb,theta0)

%% Examine the output
init = theta0.theta;
optim = sol.theta;
initResult = PStopMax(init);
% fprintf('The initial value = %d\n',)
% fprintf('The final value = %d\n',)
% Initial
GInit(1:M/2) = cos(init);
if isrow(GInit)
    GInit = GInit.';
end
GInit(M/2+1:M) = flipud(sin(init));
GInit(M+1:M+M/2) = fliplr(flipud(GInit(M/2+1:M)));
GInit(M+M/2+1:2*M) = fliplr(flipud(GInit(1:M/2)));
pInit = GInit; % The prototype filter

% Optimized
GOpt(1:M/2) = cos(optim);
if isrow(GOpt)
    GOpt = GOpt.';
end
GOpt(M/2+1:M) = flipud(sin(optim));
GOpt(M+1:M+M/2) = fliplr(flipud(GOpt(M/2+1:M)));
GOpt(M+M/2+1:2*M) = fliplr(flipud(GOpt(1:M/2)));
pOpt = GOpt; % The prototype filter

% Compute the DFT
nfft = 1024*M;
x_ax = 0:1/nfft:1-1/nfft;
x_tick = 0:1/(4*M):1-1/(4*M);
PInit = mag2db(abs(fft(pInit,nfft)));
POpt = mag2db(abs(fft(pOpt,nfft)));
figure()
plot(x_ax,PInit)
hold on
plot(x_ax,POpt,'ro')
hold off
legend('Init','Optim')
xlim([0,0.5])
xlabel('f_{norm}')
xticks(x_tick)
grid on