clc
clf
close all
clearvars

x = optimvar('x',1,2);
x0.x = [0 0];
rosenbrock = @(x)100*(x(:,2) - x(:,1).^2).^2 + (1 - x(:,1)).^2; % Vectorized function
radsqexpr = fcn2optimexpr(@disk,x);
rosenexpr = fcn2optimexpr(rosenbrock,x);
convprob = optimproblem('Objective',rosenexpr,'Constraints',radsqexpr <= 1);
[sol,fval,exitflag,output] = solve(convprob,x0)