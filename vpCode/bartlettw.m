
function s = bartlettw(u, numlag);

% BARTLETTW:
%  Computes the bartlett covariance matrix of the
%  columns of u using numlag lags.  Each lag is downweighted to
%  ensure that the resulting matrix is positive semi-definite.
%
%  This is the Newy-West covariance estimator
%
% USAGE:
%  u      - moment conditions, for example, pricing error [mR - p] (T by M)
%  numlag - number of lags considered in Newy-West estimator
%
% RESULT:
%  s - Newy-West covariance estimator
%
% © Joao Gomes, Leonid Kogan, and Lu Zhang, 2000.

[T, k]  = size(u);
mnu     = mean(u);
u       = u - ones(T, 1) * mnu;
u       = u';
s       = (u*u') /T;
iter    = zeros(k);
for lag = 1 : numlag - 1
   weight  = 1 - lag /numlag;
   iter    = (u(:,1:T-lag) * u(:,1+lag:T)' + u(:,1+lag:T)*u(:,1:T-lag)') /T;
   iter    = iter * weight;
   s       = s + iter;     
end

