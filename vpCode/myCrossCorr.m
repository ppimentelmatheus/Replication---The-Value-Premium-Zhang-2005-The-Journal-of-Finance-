
function [corrvec] = myCrossCorr(X, Y)

%-----------------------------------------------------------------------------------
% MYCROSSCORR:
%   evaluate cross correlations of two vectors X and Y
%   at -5 : 5 lags
%
% Usage:
%   X - stock return, for example
%   Y - investment return, for example
%
% OUTPUT:
%   corrvec - correlation vector containing correlations at -4 or 4 lags
%
% Modification log:
%   7/14/02: rename the file to MYCROSSCORR, different from CROSSCORR
% build-in by matlab 6.1
%-----------------------------------------------------------------------------------

corrvec = zeros(11, 1);

tmp = corrcoef(X(6 : end), Y(1 : end - 5)); % t-5
corrvec(1) = tmp(1, 2);
tmp = corrcoef(X(5 : end), Y(1 : end - 4)); % t-4
corrvec(2) = tmp(1, 2);
tmp = corrcoef(X(4 : end), Y(1 : end - 3)); % t-3
corrvec(3) = tmp(1, 2);
tmp = corrcoef(X(3 : end), Y(1 : end - 2)); % t-2
corrvec(4) = tmp(1, 2);
tmp = corrcoef(X(2 : end), Y(1 : end - 1)); % t-1
corrvec(5) = tmp(1, 2);

tmp = corrcoef(X(1 : end), Y(1 : end    )); % t
corrvec(6) = tmp(1, 2);

tmp = corrcoef(X(1 : end - 1), Y(2 : end)); % t+1
corrvec(7) = tmp(1, 2);
tmp = corrcoef(X(1 : end - 2), Y(3 : end)); % t+2
corrvec(8) = tmp(1, 2);
tmp = corrcoef(X(1 : end - 3), Y(4 : end)); % t+3
corrvec(9) = tmp(1, 2);
tmp = corrcoef(X(1 : end - 4), Y(5 : end)); % t+4
corrvec(10) = tmp(1, 2);
tmp = corrcoef(X(1 : end - 5), Y(6 : end)); % t+4
corrvec(11) = tmp(1, 2);


