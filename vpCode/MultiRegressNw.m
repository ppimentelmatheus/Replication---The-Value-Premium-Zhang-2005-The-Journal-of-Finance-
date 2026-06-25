
function [coefs, tstats, adjusted_R2, chi2, p] = MultiRegressNw(Y, X, R, q, constant)

%----------------------------------------------------------------------------
% MULTIREGRESSNW:
%   multivariate regression (Newey-West Version) including Wald test
%
% USAGE:   
%   Y --- The Independent Variable
%   X --- The Regressors excluding the Intercept
%   R --- A Matrix Providing A Linear Relationship of Coefficients
%           to Be Tested By Walt Test (p.282 in Green). R is a J by K
%           Matrix where J is the number of constraints and K is the number
%           of regressors contained in X (number of columns in X)
%   q --- The constant vector in the linear relationship
%   constant --- include constant in the regression or not
%     1 --- with intercept
%     0 --- without intercept
%
% RESULTS:
%   coefs  --- regression coefficients, intercept and slopes
%   tstats --- t statistics
%   R2     --- as itself
%   chi2   --- Walt Test statistic for linear relationship
%   p      --- p-value of the chi2 statistic
%----------------------------------------------------------------------------

[T, K]  = size(X);
[J, Jk] = size(R);

if constant == 1
   X     = [ones(T, 1) X];
end
[T, K]  = size(X);
coefs = ( inv(X'*X) )*X'*Y;

resid = Y - X * coefs;
R2 = 1 - var(resid)/var(Y);
adjusted_R2 = 1 - ((T - 1)/(T - K)) * (1 - R2);

covm = T * inv(X'*X) * bartlettw((X .* kron(resid, ones(1, K))), 12) * inv(X'*X);
stde = sqrt(diag(covm));
tstats = coefs ./stde;

if J == 0
   chi2 = [];
   p    = NaN;
elseif Jk ~= K
   disp('The dimension of the User-Provided R matrix is not conformable to that of X!')
else
   if constant == 1
      R    = [zeros(J, 1) R];
   end
   chi2 = (R*coefs - q)'* inv(R*covm*R')*(R*coefs - q);
   p    = 1 - chi2cdf(chi2, size(R, 1));
end

