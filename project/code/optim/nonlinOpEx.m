clc
clf
close all
clearvars

x = optimvar('x',1,2);
obj = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
prob = optimproblem('Objective',obj);
nlcons = x(1)^2 + x(2)^2 <= 1;
prob.Constraints.circlecons = nlcons;
showproblem(prob)
x0.x = [0 0];
[sol,fval,exitflag,output] = solve(prob,x0)