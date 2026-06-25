
function [c, t, R2, corrmtrx] = StyleTiming(Pf, Rf, Bf, Df, sx, a)

% -------------------------------------------------------------------------
% style time and explain the value spread
% -------------------------------------------------------------------------
% 
% construct regressors
[HML, SMB, ValSpr, ROEhml, mbtm, abtm] = p3by2SortEx(Pf, Rf, Bf, Df);
% 
col      = 1 : 12 : length(sx);                      % from monthly btm to annual btm

if a == 1,  % annual frequency
    HML    = TimeAggr(HML, 'm', 'y', 'f');
    ROEhml = TimeAggr(ROEhml, 'm', 'y', 'f');
    
    ValSpr = ValSpr(col);
    sx     = sx(col);
    mbtm   = mbtm(col);
    abtm   = abtm(col);
    
    % ValSpr = TimeAggr(ValSpr, 'm', 'y', 's');
    % sx     = TimeAggr(sx, 'm', 'y', 's');
    % mbtm   = TimeAggr(mbtm, 'm', 'y', 's');
    % abtm   = TimeAggr(abtm, 'm', 'y', 's');
end

xbar = mean(sx);

% style timing regressions
iota = ones(length(HML), 1);
% on constant and value spread
[c_1, t_1, R2_1, chi2, p] = MultiRegressNw(HML', [iota ValSpr'], [], [], 0); 
% on constant and earnings growth spread
[c_2, t_2, R2_2, chi2, p] = MultiRegressNw(HML', [iota ROEhml'], [], [], 0); 
% on constant and lagged earnings growth spread
[c_2lag, t_2lag, R2_2lag, chi2, p] = MultiRegressNw(HML(2:end)', [iota(1:end-1) ROEhml(1:end-1)'], [], [], 0); 
% on constant, value spread, and earnings growth spread
[c_3, t_3, R2_3, chi2, p] = MultiRegressNw(HML', [iota ValSpr' ROEhml'], [], [], 0); 
% on constant, value spread, and lagged earnings growth spread
[c_3lag, t_3lag, R2_3lag, chi2, p] = MultiRegressNw(HML(2:end)', [iota(1:end-1) ValSpr(2:end)' ROEhml(1:end-1)'], [], [], 0); 
% on constant and aggregate productivity
[c_4, t_4, R2_4, chi2, p] = MultiRegressNw(HML', [iota (sx - xbar)'], [], [], 0); 
% on constant, value spread, and market median btm
[c_5, t_5, R2_5, chi2, p] = MultiRegressNw(HML', [iota ValSpr' mbtm'], [], [], 0); 
% on constant, value spread, and market aggregate btm
[c_6, t_6, R2_6, chi2, p] = MultiRegressNw(HML', [iota ValSpr' abtm'], [], [], 0); 
% on constant, value spread, and aggregate shock
[c_7, t_7, R2_7, chi2, p] = MultiRegressNw(HML', [iota ValSpr' (sx - xbar)'], [], [], 0); 
% 
% explaining the value spread
% on constant and aggregate productivity
[c_8, t_8, R2_8, chi2, p] = MultiRegressNw(ValSpr', [iota (sx - xbar)'], [], [], 0); 
% on constant and market median btm
[c_9, t_9, R2_9, chi2, p] = MultiRegressNw(ValSpr', [iota mbtm'], [], [], 0); 
% on constant and market aggregate btm
[c_10, t_10, R2_10, chi2, p] = MultiRegressNw(ValSpr', [iota abtm'], [], [], 0); 
% on constant and lagged ROE
[c_11, t_11, R2_11, chi2, p] = MultiRegressNw(ValSpr', [iota ROEhml'], [], [], 0); 
% 
% correlation matrix
corrmtrx = corrcoef([sx(2:end)' HML(2:end)' mbtm(2:end)' ROEhml(2:end)' ROEhml(1:end-1)' ValSpr(2:end)']);

c  = [c_1; c_2; c_2lag; c_3; c_3lag; c_4; c_5; c_6; c_7; c_8; c_9; c_10; c_11];
t  = [t_1; t_2; t_2lag; t_3; t_3lag; t_4; t_5; t_6; t_7; t_8; t_9; t_10; t_11];
R2 = [R2_1; R2_2; R2_2lag; R2_3; R2_3lag; R2_4; R2_5; R2_6; R2_7; R2_8; R2_9; R2_10; R2_11];