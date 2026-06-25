
function [LV1fac, LV2fac, LV3fac] = liewVass(simGDPg, simHML, simSMB, simRm)

% LIEWVASS:
%   Conducting on the simulation data Liew and Vassalou (2000) time series regressions 
% (all regressions are at quarterly frenquency):
% 
% 1. Three univariate regressions: GDP growth on market, HML, and SMB
%     GDP growth(t, t+4) = a + b factor return(t-4, t) + residual
% 
% 2. Two bivariate regressions: GDP growth on market and SMB or on market and HML
%     GDP growth(t, t+4) = a + b market(t-4, t) + c factor return(t-4, t) + residual
% 
% 3. One three-factor regressions:
%     GDP growth(t, t+4) = a + b market(t-4, t) + c HML(t-4, t) + d SMB(t-4, t) + residual
% 
% USAGE: 
%   simGDPg -- 1 by T monthly GDP growth, where T is a multiply of 360
%   simHML  -- 1 by T monthly HML return, where T is a multiply of 360
%   simSMB  -- 1 by T monthly SMB return, where T is a multiply of 360
%   simRm   -- 1 by T monthly mkt return, where T is a multiply of 360
% where the timing is exactly that simGDPg leads the other three factor returns by 12 months
% 
% OUTPUT:
%   LV1fac  -- 3 by 5 matrix (on each row [a b ta tb r2]; on each column [market SMB HML]')
%   LV2fac  -- 2 by 7 matrix (on each row [a b c ta tb tc r2]; on each column [SMB HML]')
%   LV3fac  -- 1 by 9 matrix ([a b c d ta tb tc td r2])
% where ta/b/c/d denotes t-stats and r2 denotes goodness of fit coefficient
%
% © Lu Zhang, Inc. 2001


% transform simulated data to nsim samples
T    = length(simGDPg);
nsim = T/360;           
GDPg = reshape(simGDPg - 1, 360, nsim)';
HML  = reshape(simHML,  360, nsim)';
SMB  = reshape(simSMB,  360, nsim)';
Rm   = reshape(simRm,   360, nsim)';

% cut the sample size to 20 years worth each, 
% Liew and Vassalou (2000) use 19 years worth of sample (the extra one year used be used to construct 
%    overlapping annualized returns in quarterly frenquency)
% note sample at this point is still in monthly frequency
GDPg = GDPg(:, 1:240);
HML  =  HML(:, 1:240);
SMB  =  SMB(:, 1:240);
Rm   =   Rm(:, 1:240);

% now transform the sample to annualized values in quarterly frequency
% 
% step 1: transform monthly sample to non-overlapping quarterly sample (nsim by 80 matrix)
GDPg = squeeze(sum(reshape(GDPg', [3, 80, nsim]), 1))';
HML  = squeeze(sum(reshape(HML', [3, 80, nsim]), 1))';
SMB  = squeeze(sum(reshape(SMB', [3, 80, nsim]), 1))';
Rm   = squeeze(sum(reshape(Rm', [3, 80, nsim]), 1))';
% 
% step 2: now construct overlapping quarterly sample in annual terms --- the last year of data will be used up here
%    results are nsim by 76 matrix (76/4 = 19 years)
samGDPg = zeros(nsim, 76);
samHML  = zeros(nsim, 76);
samSMB  = zeros(nsim, 76);
samRm   = zeros(nsim, 76);
% 
for n = 1 : nsim        % the last period is redundant
    samGDPg(n, :) = GDPg(n, 1:76) + GDPg(n, 2:77) + GDPg(n, 3:78) + GDPg(n, 4:79);
     samHML(n, :) =  HML(n, 1:76) +  HML(n, 2:77) +  HML(n, 3:78) +  HML(n, 4:79);
     samSMB(n, :) =  SMB(n, 1:76) +  SMB(n, 2:77) +  SMB(n, 3:78) +  SMB(n, 4:79);
      samRm(n, :) =   Rm(n, 1:76) +   Rm(n, 2:77) +   Rm(n, 3:78) +   Rm(n, 4:79);
end

% Univariate regressions
LV1fac = zeros(3, 5, nsim);
for n = 1 : nsim
    % univariate regression of GDP growth on market
    LV1fac(1, :, n) = RegressNw(samGDPg(n, :)',  samRm(n, :)');
    % univariate regression of GDP growth on SMB 
    LV1fac(2, :, n) = RegressNw(samGDPg(n, :)', samSMB(n, :)');
    % univariate regression of GDP growth on HML 
    LV1fac(3, :, n) = RegressNw(samGDPg(n, :)', samHML(n, :)');
end
% 
% Bivariate regressions
LV2fac = zeros(2, 7, nsim);
for n = 1 : nsim
    % GDP growth on market and SMB
    LV2fac(1, :, n) = RegressNw(samGDPg(n, :)', [samRm(n, :)' samSMB(n, :)']);
    % GDP growth on market and HML
    LV2fac(2, :, n) = RegressNw(samGDPg(n, :)', [samRm(n, :)' samHML(n, :)']);
end
% 
% 3 factor regression
LV3fac = zeros(1, 9, nsim);
for n = 1 : nsim
    LV3fac(1, :, n) = RegressNw(samGDPg(n, :)', [samRm(n, :)' samHML(n, :)' samSMB(n, :)']);
end

% average over samples to obtain final results
LV1fac = squeeze(mean(LV1fac, 3));
LV2fac = squeeze(mean(LV2fac, 3));
LV3fac = mean(LV3fac, 3);
    
    
% ------------------------------ subfunction declaration -------------------------
% 
function [results] = RegressNw(Y, x)

% RegressNw:%   Linear regression with Newey and West adjusted t-stats. 
% Y is the dependent variable and x is the independent variable: both are column vectors
% 
% OUTPUT: 
%   coef   -- size(x, 2) + 1 row vector of regression coefficients
%   tstats -- size(x, 2) + 1 row vector of t-stats for intercept adjusted by Newey and West procedure%   r2adj  -- adjusted r2
% initializing
T    = length(Y);
X    = [ones(length(x), 1) x];
% coefficients
coef = inv(X'*X)*X'*Y;
residual = Y - X*coef;
coef = coef';
% 
% tstats
if 0
covm     = T*inv(X'*X)*bartlettw(residual, 3)*inv(X'*X);tstats   = sqrt(T)*coef ./sqrt(diag(covm)');
else
covm     = (residual'*residual /(T - size(X, 2))) * inv(X'*X);
tstats   = coef ./sqrt(diag(covm)');
end    
% 
% adjusted r2
R2       = 1 - var(residual)/var(Y);
r2adj    = 1 - (T - 1)/(T - 2) *(1 - R2);

results  = [coef tstats r2adj];  

% ------------------------------ subfunction declaration -------------------------
% 
function [s] = bartlettw(u, numlag);

% BARTLETTW:
%  Computes the bartlett covariance matrix of the
%  columns of u using numlag lags.  Each lag is downweighted to
%  ensure that the resulting matrix is positive semi-definite.
%
%  This is the Newy-West covariance estimator
%
% USAGE:
%  u - moment conditions, for example, pricing error [mR - p] (T by M)
%  numlag - number of lags considered in Newy-West estimator
%
% RESULT:
%  s - Newy-West covariance estimator
%
% © Joao Gomes, Amir Yaron, and Lu Zhang, 2000.

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

