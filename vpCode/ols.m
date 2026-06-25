 
function [b, tstats, R2, sighat, resid] = ols(Y, X)

% OLS
%   Simple linear regression routine used in predicting future aggregate state or variables

b      = inv(X'*X)*X'*Y;
resid  = Y - X*b;
sighat = std(resid);

R2     = 1 - var(resid)/var(Y);
covm   = (resid'*resid /(length(Y) - size(X, 2))) * inv(X'*X);
tstats = b ./sqrt(diag(covm));

