% PANIE:
%   Simulate Panel Data from Cross-Sectional Firm Distribution and Perform Fama-French (1992, 1993)
% and Exact Regression and Portfolio Grouping as well as HML and SMB summary statistics
% 
% See SSIE for simulating the stationary distribution
%
% © Lu Zhang, Inc. 2001

clear; clc; format short; format compact; warning off

load Params                                                   % parameters
load distrSS                                                  % stationary distribution of firms
load vfi3Mat 
load coefIS

% reshape to facilatate interpolation
optK = reshape(optK, [nk nh nx nz]);
V0   = reshape(V0, [nk nh nx nz]);

% initializing
N        = 5000;                                              % number of firms in the sample
Ts       =  421;                                              % number of periods (month)
nsim     =   20;                                              % number of repetitive simulations
% tbl99  = zeros(3, 3, nsim);
tbl99    = zeros(3, 9, nsim);
ffMtrx   = zeros(2, 9, nsim);
corrMtrx = zeros(1, 3, nsim);
table1   = zeros(11, 11, nsim);  table2 = zeros( 4, 12, nsim);  table3 = zeros( 4, 12, nsim);
table4   = zeros( 4, 12, nsim);  tbl93  = zeros(30, 10, nsim);
% some intermediate variables concerning asymmetric risk story
simSMB   = [];    simHML   = [];    simSH    = [];    simBH    = [];    simSL    = [];    simBL    = [];
simEP    = [];    simRm    = []; 
simx     = [];    simrf    = [];    simBTM   = [];  
simCORR  = [];    simGDPg  = [];    simRmLV  = [];    simRmSH  = [];    simSDF   = [];
corrTimeSeries = zeros(nsim, 1);
% shares
iyr   = zeros(nsim, 1); 
ikr   = zeros(nsim, 1);
theta = zeros(nsim, 1);
dyr   = zeros(nsim, 1);
fyr   = zeros(nsim, 1);
negI  = zeros(nsim, 1);
posI  = zeros(nsim, 1);

% start repetitive simulation consecutively along time dimension
tic;
for is   = 1 : nsim
%--------------------------------------------------

    % aggregate shocks 
    [sx, xinlnew] = CspSimu(xinl, xbar, rhox, stdx, Ts);
    xinl  = xinlnew;
    % 
    % stochastic discount factor
    SDF   = beta*exp((gamA + gamB*(sx(1:end-1) - xbar)).*(sx(1:end-1) - sx(2:end)));
    
    % real interest rate
    srf   = interp1(x, rf, sx, 'linear', 'extrap');    srf(end) = [];

    % simulation MEX routine
    [Pf, Bf, Df, Rf, In, iyr(is), ikr(is), theta(is), dyr(is), fyr(is), Rm, GDPg, kd0new, zd0new] = ...
       panIEfcn(kd0, zd0, optK, V0, k, sx, h, nh, x, nx, z, N, Ts, alpha, alp1, alp2, alp3, delta, eta, f, gP, gN, istar, rhoz, stdz);
    Pf(:, end) = [];     
    Bf(:, end) = [];    
    Df(:, end) = []; 
    In(:, end) = [];
    negInd     = find(In < 0);
    posInd     = find(In > 0);
    negI(is)   = mean(In(negInd) ./Bf(negInd));
    posI(is)   = mean(In(posInd) ./Bf(posInd));
    clear In negInd posInd
    
    % initial distribution for next iteration
    kd0 = kd0new;
    zd0 = zd0new;
    % 
    % quality control
    if min(min(Pf)) < 0
        fprintf('negative firm value!\n');
        Pf = max(Pf, 1e-8);
    end
    
    % perform fama and french (1992) cross-sectional regression tests 
    [ffMtrx(:, :, is), corrMtrx(:, :, is), Portmtrx, postbeta, FirmIndex] = FF92(Pf, Rf - 1, Bf, Rm - 1, N);
    % 
    % Perform fama-french portfolio grouping to report summary statistics
    [table1(:, :, is), table2(:, :, is), table3(:, :, is), table4(:, :, is)] = ...
                     shortPort(Rf - 1, Pf, Bf, Rm - 1, Portmtrx, postbeta, FirmIndex, N);
    % 
    % House Keeping
    clear Portmtrx postbeta FirmIndex 
    % 
    % summary statistics of HML and SMB
    [tbl99(:, :, is), SMB, HML, SH, BH, SL, BL] = ValPrem(Pf, Rf, Bf, Rm, srf);
    % 
    % perform fama and french (1993) time-series analysis
    tbl93(:, :, is) = FF93(Pf, Rf, Bf, Rm, srf, SMB, HML);
    % 
    % perform fama and french (1995) data crunching
    [roevec(:, :, is), btmvec(:, :, is), ROElow(:, :, is), ROEhigh(:, :, is)] = FF95(Pf, Bf, Df);

    % save data for time series regressions
    index   = find(isnan(HML));
    if length(index) == 0    
        simGDPg = [simGDPg GDPg(61 : end)];   % save GDP growth with timing exactly one year leading factor returns
        simSMB  = [simSMB SMB];    simHML  = [simHML HML];
        simSH   = [simSH SH];      simBH   = [simBH BH];    simSL   = [simSL SL];      simBL   = [simBL BL];
        simEP   = [simEP  log(Rm) - log(srf)];
        simRm   = [simRm Rm];
        simRmLV = [simRmLV log(Rm(61 : end)) - log(srf(61 : end))];
        simRmSH = [simRmSH Rm(61 : end)];
        simx    = [simx    sx(61 : end - 1)];
        simSDF  = [simSDF  SDF(61 : end)];
        simrf   = [simrf  srf];
        simBTM  = [simBTM sum(Bf)./sum(Pf)];

        % cross-sectional and time-series average correlation between btm and returns
        tmp = corrcoef(simBTM, simRm); 
        corrTimeSeries(is) = tmp(1, 2);
        btm  = Bf./Pf;
        btm  = [zeros(N, 6)    btm(:, 1 : 420 - 6)];
        corrbtm  = zeros(1, 360);
        corrsize = zeros(1, 360);
        for t = 61 : 420
            tmp              = corrcoef(btm(:, t), Rf(:, t));
            corrbtm(t - 60)  = tmp(1, 2);
            tmp              = corrcoef(Pf(:, t), Rf(:, t));
            corrsize(t - 60) = tmp(1, 2);
        end
        simCORR = [simCORR [mean(corrbtm); mean(corrsize)]];
    else
        fprintf('NaN observed...\n');
    end
    
    % visual
    if mod(is, 50) == 0,
        fprintf('simulation %2.0f is done within %6.2f seconds.\n', [is toc]);
    end
    
%-----------------------------------------------------------
end  % nsim

% NOW REPORT EMPIRICAL RESULTS

% PART 1: TIME SERIES STATISTICS
% 
% equity premium
meanEP = mean(simEP)*12
stdEP  = std(simEP)*sqrt(12)
% 
% real interest rate
meanrf = mean(simrf)^12
stdrf  = std(simrf)*sqrt(12)
% 
% btm
mBTM   = mean(simBTM)
sBTM   = std(simBTM)
% 
% aggregate ratios
shares = mean([iyr ikr theta dyr fyr negI posI])

% PART 2: TIME SERIES REGRESSIONS
% 
% Pontiff and Schall (1999)
btmPanel = ponSchall(simRm, simBTM)
% 
% on ff 93 time-series regression
tbl93    = mean(tbl93, 3)
% 
% Lewellen (1999): Also want to show value premium measured in this alternative way is higher than Gomes, Kogan, and Zhang (2001)
%   To be inserted ...
% 
% Liew and Vassalou (2000)
[LV1fac, LV2fac, LV3fac] = liewVass(simGDPg, simHML, simSMB, simRmLV)
% 
corrHmlMkt = corrcoef(simHML', simRmSH')

% PART 3: THE CROSS-SECTION
% 
% Fama and French (1992) portfolio grouping
table1  = mean(table1, 3)
table2  = mean(table2, 3)
table3  = mean(table3, 3)
table4  = mean(table4, 3)
% 
% Summary statistics of unconditional size and value premium 
tbl99    = mean(tbl99, 3)   
sizePrem = [mean(simSMB)*12*100 std(simSMB)*100*sqrt(12)]
valPrem  = [mean(simHML)*12*100 std(simHML)*100*sqrt(12)]
% 
% on cross-sectional regression
ff      = mean(ffMtrx, 3)
mcorr   = mean(corrMtrx, 3)
corrcoef([simx' simHML' simSMB' simSH' simBH' simSL' simBL'])
% 
% cross-sectional average of time series correlation between btm and return
corrTS = mean(corrTimeSeries)
corr   = mean(simCORR, 2)
figure(10); plot(simCORR(1,:)); xlabel('simulation'); ylabel('corr');
hold on; plot(simCORR(2,:), 'r'); xlabel('simulation'); 
legend('corr(btm, Rf)', 'corr(size, Rf)');
% It seems that both the signs and magnitudes of the two correlation series are sensible

% FF 95
roevec  = mean(roevec, 3);
btmvec  = mean(btmvec, 3);
ROElow  = mean(ROElow, 3);
ROEhigh = mean(ROEhigh, 3);
[Junk]  = plotFF95(roevec, btmvec, ROElow, ROEhigh); 

% PART 4: Lettau and Ludvigson (2001) table on asymmetric risk
tableLL = zeros(5, 3);
good    = find(simx > xbar + (1/2)*stdx/sqrt(1 - rhox^2));
bad     = find(simx < xbar - (2)*stdx/sqrt(1 - rhox^2));
% ALL STATES
xseries = simSDF;   % regress against SDF
% xseries = simx;
XX = [ones(length(xseries), 1) xseries'];
[b, tstats, R2] = ols(simSH', XX);
tableLL(1, 1)   = b(2); 
[b, tstats, R2] = ols(simBH', XX);
tableLL(2, 1)   = b(2); 
[b, tstats, R2] = ols(simSL', XX);
tableLL(3, 1)   = b(2); 
[b, tstats, R2] = ols(simBL', XX);
tableLL(4, 1)   = b(2); 
[b, tstats, R2] = ols(simHML', XX);
tableLL(5, 1)   = b(2); 
% GOOD STATES
XX = [ones(length(xseries(good)), 1) xseries(good)'];
[b, tstats, R2] = ols(simSH(good)', XX);
tableLL(1, 2)   = b(2); 
[b, tstats, R2] = ols(simBH(good)', XX);
tableLL(2, 2)   = b(2); 
[b, tstats, R2] = ols(simSL(good)', XX);
tableLL(3, 2)   = b(2); 
[b, tstats, R2] = ols(simBL(good)', XX);
tableLL(4, 2)   = b(2); 
[b, tstats, R2] = ols(simHML(good)', XX);
tableLL(5, 2)   = b(2); 
% BAD STATES
XX =  [ones(length(xseries(bad)), 1) xseries(bad)'];
[b, tstats, R2] = ols(simSH(bad)', XX);
tableLL(1, 3)   = b(2); 
[b, tstats, R2] = ols(simBH(bad)', XX);
tableLL(2, 3)   = b(2); 
[b, tstats, R2] = ols(simSL(bad)', XX);
tableLL(3, 3)   = b(2); 
[b, tstats, R2] = ols(simBL(bad)', XX);
tableLL(4, 3)   = b(2); 
[b, tstats, R2] = ols(simHML(bad)', XX);
tableLL(5, 3)   = b(2); 

tableLL  % = abs(tableLL)