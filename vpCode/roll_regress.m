
function [beta, R2] = roll_regress(Rp, Rm, len)

%------------------------------------------------------------------------
% ROLL_REGRESS: Estimate conditional beta by rolling window estimation
% 
% INPUT:
%   Rp  - portfolio EXCESS returns: T by 1
%   Rm  - market EXCESS return: T by 1
%   len - length of the window used in rolling regression 
% 
% OUTPUT:
%   beta - time series of conditional beta: T - len
%   R2   - time series of R2: T - len
% 
% NOTE: Rm can also be other business cycle variables
%------------------------------------------------------------------------

len_p = length(Rp);
len_m = length(Rm);

if len_p ~= len_m, error('The length of Rp and Rm are not equal!'); end
if len_p < len, error('Window length is too big or Rp and Rm time series is too short!'); end

beta = zeros(len_p - len, 1);
R2   = zeros(len_p - len, 1);
for t = len + 1 : len_p
    % use current and lagged-one-period market returns
    % X = [ones(len, 1) Rm(t - len + 1 : t) Rm(t - len : t - 1)]; 
    % tmp = inv(X'*X)*X'*Rp(t - len + 1 : t);
    % beta(t - len) = sum(tmp(2 : 3));
    % R2(t - len) = var(X*tmp)/var(Rp(t - len + 1 : t));
    % use only current market returns
    X = [ones(len, 1) Rm(t - len + 1 : t)]; 
    tmp = inv(X'*X)*X'*Rp(t - len + 1 : t);
    beta(t - len) = tmp(2);
    R2(t - len) = var(X*tmp)/var(Rp(t - len + 1 : t));
end
    