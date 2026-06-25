
function [em, sigm, Sp, lambda] = lambda_m(x, beta, gamA, gamB, xbar, rhox, stdx)

% LAMBDA_M:
%    Annualized Sharpe ratio as a function of aggregate shock given other parameter values where pricing kernel:
% 
%                 log(M_{t,t+1}) = log(beta) + (gamA + gamB*(x_t - xbar))*(x_t - x_{t+1})
% 
%   em   - E_t[M_{t+1}]
%   sigm - \sigma_t[M_{t+1}]
%   Sp   - sigm/em
%   lambda - sigm^2/em
% 
% © Lu Zhang, Inc. 2001.

mu  = (gamA + gamB*(x - xbar))*(1 - rhox).*(x - xbar);
sig = stdx*(gamA + gamB*(x - xbar));

em   = beta*exp(mu + (1/2)*sig.^2);
sigm = beta*exp(mu).*sqrt( exp(sig.^2) .* (exp(sig.^2) - 1) );
Sp   = sigm./em;
lambda = sigm.^2./em;

% etmp   = exp((theta_0 + theta_1*(x - xbar)).^2 * stdx^2);
% varm   = etmp.*(etmp - 1); 
% em     = exp((1/2)*(theta_0 + theta_1*(x - xbar)).^2 * stdx^2); 
% lambda = varm./em; 