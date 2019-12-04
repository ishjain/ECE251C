% Optimize the Kaiser window such that the power complementary condition is
% satisfied to the maximum extent
clc
clf
close all
clearvars

%% Prototype filter parameter definition
M = 16;
m = 1;
N = 2*m*M-1;
L = N+1;

%% Brute-force search
M = 16;
m = 1;
diff = @(alpha) diffFromPowComp(alpha,M,m);
alpha = 22:0.01:60;
for i = 1:length(alpha)
    result(i) = diff(alpha(i));
end

[minResult,minidx] = min(result);
alpha(minidx)

%% Optimization definition
% % Optimization variable
% alpha = optimvar('alpha');
% alpha0.alpha = 40; % Initial value
% % Objective function to optimize
% % M = 16;
% % m = 1;
% % diff = @(alpha) diffFromPowComp(alpha,M,m);
% objExpr = fcn2optimexpr(@diffFromPowComp,alpha);
% % Constraints
% minVal = @(alpha) min(alpha);
% minConstr = fcn2optimexpr(minVal,alpha);
% % Define the optim problem
% convProb = optimproblem('Objective',objExpr,'Constraints',minConstr>=21);
% [sol,fval,exitflag,output] = solve(convProb,theta0)

%% Check the output
% fprintf('The initial value = %d\n', diff(alpha0.alpha));
% fprintf('The final value = %d\n', fval);