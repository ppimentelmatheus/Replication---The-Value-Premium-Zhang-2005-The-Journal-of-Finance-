
function ySharpe = ySharpeFcn(x, theta_0, theta_1, xbar, stdx)

% YSHARPEFCN:
%    Annualized Sharpe ratio as a function of aggregate shock given other parameter values where pricing kernel:
% 
%                 log(M_{t,t+1}) = log(beta) + (theta + theta_1*(x_t - xbar))*(x_t - x_{t+1})
% 
% © Lu Zhang, Inc. 2001.

etmp    = exp((theta_0 + theta_1*(x - xbar)).^2 * stdx^2);
ySharpe = sqrt(12)*sqrt(etmp.*(etmp - 1)) ./exp((1/2) * (theta_0 + theta_1*(x - xbar)).^2 * stdx^2); 
