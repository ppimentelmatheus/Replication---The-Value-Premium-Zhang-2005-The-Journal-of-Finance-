% SSIE:
%   Simulate the stationary distribution of firms until stationary distribution
%
% OUTPUT: (Saved in distrSS.mat)
%   kd   - Cross-Sectional Capital Stock
%   zd   - Cross-Sectional Idiosyncratic Shock
%
% See SIMPAN for simulating cross-sectional statistics from stationary distribution
% 
% © Lu Zhang, Inc. 2001

clear all; clc; format short; format compact;

load Params
load vfi3Mat

% reshape optK to facilitate linear interpolation
optK   = reshape(optK, [nk nh nx nz]);

% Initialize the economy
N      =  5000;                                                 % number of firms in the sample
Ts     = 10000;                                                 % number of periods (month) to achieve stationarity

if 0
    load distrSS
    % aggregate shock from continuous state space
    [sx, xinl] = CspSimu(xinl, xbar, rhox, stdx, Ts);
else
    % initializing the firm distribution ('d' denotes values in sampling distribution)
    kd0    = ones(N, 1);                                        % cross-sectional current period capital stock
    % idiosyncratic shock from continuous state space
    shock  = randn(N/2, 1);                                     % impose law of large numbers using antithetic methods
    zd0    = max(-3.5*(stdz/sqrt(1 - rhoz^2)), min(3.5*(stdz/sqrt(1 - rhoz^2)), (stdz/sqrt(1 - rhoz^2)) * [shock; -shock]));
    % aggregate shock from continuous state space
    [sx, xinl] = CspSimu(xbar, xbar, rhox, stdx, Ts);
end
xinl
% call the bottle-neck simulation MEX routine 
[kpd, zpd, simh] = ssIEfcn(kd0, zd0, sx, optK, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz);
kd0  =  kpd;              % starting value --- k(1) in PANIE_P routine
zd0  =  zpd;              % starting value --- z(1) in PANIE_P routine
% sh0  = simh(end);       % starting value --- simh(0) in PANIE_P routine 
% sx0  = sx(end);         % starting value --- sx(0) in PANIE_P routine --- note sx(1) is saved as xinl
% 
% Not necessary to keep track of sh0 and sx0 since log output price can be calculated easily
% from simulated firm matrix; high R2p only validates the use of approximate law of motion in
% solving firm's problem. 

% compute some simple statistics on stationary distribution without simulating panels
stats = [ mean(kpd) std(kpd) skewness(kpd) kurtosis(kpd);
          mean(zpd) std(zpd) skewness(zpd) kurtosis(zpd) ] 

% figure(51); hist(kd0, 50)
% figure(52); plot(simh)


% save stationary distribution as the starting point of SIMPAN routine
save distrSS xinl kd0 zd0 

panIE900