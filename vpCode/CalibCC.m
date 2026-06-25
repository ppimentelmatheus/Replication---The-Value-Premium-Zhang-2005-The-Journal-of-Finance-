% ----------------------------------------------------
% CALIBCC:
%   Calibrate the economy in the monthly frequency
% 
% © Lu Zhang, Inc. 2001.
% ----------------------------------------------------

clear all; format short; format compact; clc

% pricing kernel: log(M_{t,t+1}) = log(beta) + (gamA + gamB*(x_t - xbar))*(x_t - x_{t+1})
beta    =   .994;                                    % .994 AS THE BENCHMARK WITH GAMMA1 = -1000
gamA    =     50;
gamB    =  -1000;                                    % -1,000 as benchmark

% aggregate demand: the inverse of price elasticity of demand (in [0 1] --- eta = 0 reduces to P = 1)
eta     = .50;

% technology parameter values
alpha   = .30;			                             % Capital's Share in Ouput
delta   = .01;			                             % Depreciation Rate of Capital: delta  = 12\% per annum

% parameters in asymmetric adjustment cost specification
gP      =  15;                                       % benchmark: gP = 15
gN      = 150;                                       % benchmark: gN = 150
istar   =   0;
cfrac   = (gP/2)*(delta - istar)^2/delta             % steady state share of adjustment cost in investment

% idiosyncratic shocks -- Rouwenhorst (1995) method
nz      = 15;
zbar    = 0;										 % mean
rhoz    = 0.97;                                      % 
stdz    = 0.10;                                      % benchmark level is .10
zdev    = 2*stdz/sqrt((1 - rhoz^2)*(nz - 1));        % do ''help rouwTrans'' to understand the formular
[Qz, z] = rouwTrans(rhoz, zbar, zdev, nz);
zmin    = min(z);
zmax    = max(z);
imz     = ceil(nz/2);								 % index of the average shock
uQz     = Qz^1000;                                   % unconditional distribution of z
uQz     = uQz(:, 1);

% simulation parameters
N       =  5000;                                     % number of firms 
Ts      = 11000;                                     % number of periods in Krusell and Smith simulation

% aggregate shocks -- Rouwenhorst (1995) method
nx      =  11;    								     % number of grid points
rhox    = .95^(1/3);                                 % calibrated to match rhox = .95  quarterly
stdx    = .007/3;                                    % calibrated to match stdx = .007 quarterly
% Calibrating xbar to normalize ks = 1: Adjusting xbar for Jensen's Inequality term associated with z shock
xbar    = (1/(1 - eta))*log((1 + gP*(delta - istar) - beta*((gP/2)*(delta - istar)*(delta + istar) + ...
    (1 - delta)*(1 + gP*(delta - istar))))/(alpha*beta*exp((1/2)*(1 - eta)^2*stdx^2/(1 - rhox^2) + ...
    (1/2)*(1 - eta)^2*stdz^2/(1 - rhoz^2))))
xdev    = 2*stdx/sqrt((1 - rhox^2)*(nx - 1));    
[Qx, x] = rouwTrans(rhox, xbar, xdev, nx);
xmax    = max(x); 
xmin    = min(x);
imx     = ceil(nx/2);
volr    = stdz/stdx

% combined state space and transition matrix for exogenous shocks: x changes faster than z
Qzx     = kron(Qz, Qx);

% average Sharpe ratio
mSharpe = sqrt(exp(gamA^2*stdx^2) * (exp(gamA^2*stdx^2) - 1)) /exp(1/2*gamA^2*stdx^2);
ySharpe = mSharpe * sqrt(12)

% Sharpe ratio across business cycle
% fplot('ySharpeFcn', [xmin xmax], 1e-3, [], [], gamA, gamB, xbar, stdx); grid on; hold on;
% fplot('ySharpeFcn', [xmin xmax], 1e-3, [], [], gamA, 0, xbar, stdx);
if 0
plotSp  = ySharpeFcn(x, gamA, gamB, xbar, stdx);
plot(x, plotSp, 'LineWidth', 2); grid on;
set(gca, 'FontS', 15);
xlabel('Aggregate Shock', 'FontS', 15); ylabel('Market Price of Risk', 'FontS', 15)
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\Sharpe.eps
end

% fixed cost of production
f    = .0365                      % benchmark 0.0365

% some steady-state relations
% ys   = exp((1 - eta)*xbar + (1/2)*stdx^2/(1 - rhox^2) + (1/2)*stdz^2/(1 - rhoz^2))  % - eta*alpha*log(N) + (1/2)*stdz^2/(1 - rhoz^2)          % revenue
% is   = delta;                                                                       % investment
% ds   = ys - f - is - (gP/2)*(delta - istar)^2;                                      % dividend
% aggr = ds + is + f + (gP/2)*(delta - istar)^2 - ys;                                 % check aggregation
% iyr  = is/ys                                                                        % investment/output 
% dyr  = ds/ys                                                                        % dividend/output
% ayr  = (gP/2)*(delta - istar)^2/ys
% fyr  = f/ys

% real interest rate solved on the grid of x
rf   = getRfcc(x, Qx, beta, gamA, gamB)';
if 1
% aggregate shocks on discrete state space
[sx, Junk]   = CspSimu(xbar, xbar, rhox, stdx, 205000);
sx(1 : 5000) = [];
% real interest rate
rfs = interp1(x, rf, sx, 'linear', 'extrap');
% rfs = interp1(x, rf, sx, 'nearest');
mrf   = mean(log(rfs))*12
vrf   = std(rfs)*sqrt(12)
% clear sx rfs mrf vrf
end

% Construct grid for capital stock
kmin  = 0.01;
kmax  =   10;

if 0
    % construct grid for capital stock --- use equal-spaced grid
    nk = 50;
    k  = linspace(kmin, kmax, nk)';
else
    % construct grid for capital stock --- use nonequal-spaced grid via McGrattan method
    next  = 1;
    ik    = 1;
    k     = kmin;
    % recursive construction: some experiment shows that 50 grid points are enough
    while next < kmax
        ik    = ik + 1; 
        % next = k(ik - 1) + .005*exp(0.07475*(ik - 2));                      %  50 Grid Points (kmax = 2.5)
        next = k(ik - 1) + .005*exp(0.28165*(ik - 2));                        %  25 Grid Points (kmax =  10)
        % next = k(ik - 1) + .005*exp(0.11445*(ik - 2));                      %  50 Grid Points (kmax =  10)
        % next = k(ik - 1) + .005*exp(0.302875*(ik - 2));                     %  25 Grid Points (kmax =  15)
        % next = k(ik - 1) + .005*exp(0.12475*(ik - 2));                      %  50 Grid Points (kmax =  15)
        % next = k(ik - 1) + .005*exp(0.32935*(ik - 2));                      %  25 Grid Points (kmax =  25)
        % next = k(ik - 1) + .005*exp(0.137575*(ik - 2));                     %  50 Grid Points (kmax =  25)
        % next = k(ik - 1) + .005*exp(0.1546275*(ik - 2));                    %  50 Grid Points (kmax =  50)
        % next = k(ik - 1) + .005*exp(0.0667*(ik - 2));                       % 100 Grid Points (kmax =  50)
        if next < kmax, k(ik) = next; end
    end
    k       = k';
    if 1 % put ks = 1 into the grid by brutal force
        [Junk, index] = min(abs(k - 1));   clear Junk
        if k(index) < 1
            k   = [k(1:index); 1; k(index+1:end)];
            imk = index + 1;
        elseif k(index) > 1
            k   = [k(1:index-1); 1; k(index:end)];
            imk = index;
        end
    end
end
clear next ik 
kmax  = max(k)
nk    = length(k)

% construct grid of next period capital
nkp   = 5000
kp    = linspace(kmin, kmax, nkp)';

% construct grid of log output price --- the end points are obtained via simulation
nh   = 5;
hmin = 2.75;
hmax = 3.25;
h    = linspace(hmin, hmax, nh)';
imh  = ceil(nh/2);

% save parameters for later use
save Params