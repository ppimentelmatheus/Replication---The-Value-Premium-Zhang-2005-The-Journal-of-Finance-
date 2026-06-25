% PANIE900:
%   Simulate Panel Data from the Cross-Sectional Distribution of Firms and construct Fama-French 25
% portfolio, 10 B/M portfolio, and 10 Size portfolio
% 
% See SSIE for simulating the stationary distribution
%
% © Lu Zhang, Inc. 2001

clear; clc; format short; format compact; warning off

load ParamsR3                                                 % parameters
load distrSSR3                                                % stationary distribution of firms
load vfi3MatR3
load coefISR3

% reshape to facilatate interpolation
optK = reshape(optK, [nk nh nx nz]);
V0   = reshape(V0, [nk nh nx nz]);

% initializing
N        = 5000;                                              % number of firms in the sample
Ts       =  900 + 1;                                          % number of periods (month)
nsim     =   10;                                              % number of repetitive simulations

% shares
iyr   = zeros(nsim, 1); 
ikr   = zeros(nsim, 1);
theta = zeros(nsim, 1);
dyr   = zeros(nsim, 1);
fyr   = zeros(nsim, 1);

negI  = zeros(nsim, 1);
posI  = zeros(nsim, 1);

btmPanel  = zeros(2, 5, nsim);
rovPanel  = zeros(2, 5, nsim);
equiPanel = zeros(4, 7, nsim);

% sHML  = []; 
% sMKT  = [];
% simrf = [];
% sbtm  = [];

if 1 % long initialization for moments of sorting portfolios
mbtm      = zeros(10, nsim);
stdbtm    = zeros(10, nsim);    
msize     = zeros(10, nsim);
stdsize   = zeros(10, nsim);
mhml      = zeros(nsim, 1);
vhml      = zeros(nsim, 1);
msmb      = zeros(nsim, 1);
vsmb      = zeros(nsim, 1);
alp_size  = zeros(10, nsim);
talp_size = zeros(10, nsim);
bet_size  = zeros(10, nsim);
tbet_size = zeros(10, nsim);
alp_btm   = zeros(10, nsim);
talp_btm  = zeros(10, nsim);
bet_btm   = zeros(10, nsim);
tbet_btm  = zeros(10, nsim);
alp_hml   = zeros(nsim, 1);
talp_hml  = zeros(nsim, 1);
bet_hml   = zeros(nsim, 1);
tbet_hml  = zeros(nsim, 1);
alp_smb   = zeros(nsim, 1);
talp_smb  = zeros(nsim, 1);
bet_smb   = zeros(nsim, 1);
tbet_smb  = zeros(nsim, 1);
end

% percbtm   = zeros(10, nsim);
% percsize  = zeros(10, nsim);
% perchml   = zeros(nsim, 1);
% percsmb   = zeros(nsim, 1);

if 0 % beta-premium analysis
out_btm = zeros(11, 4, nsim);
end

if 1
% style timing initialization
c_m        = zeros(31, nsim);     % monthly frequency
t_m        = zeros(31, nsim);
R2_m       = zeros(13, nsim);
corrmtrx_m = zeros(6, 6, nsim);
c_a        = zeros(31, nsim);     % annual frequency
t_a        = zeros(31, nsim);
R2_a       = zeros(13, nsim);
corrmtrx_a = zeros(6, 6, nsim);
end

if 0
% ff95 initialization
roevec  = zeros(11, 2, nsim);
btmvec  = zeros(11, 2, nsim);
ROElow  = zeros(nsim, 75);
ROEhigh = zeros(nsim, 75);
end

if 0
% initialization of book portfolio statistics
mvec   = zeros(10, nsim);
bvec   = zeros(10, nsim);
sigvec = zeros(10, nsim);
end

% start repetitive simulation consecutively along time dimension
tic;

if 0
% average annualized individial volatility
avevol_1 = zeros(nsim, 1);
% average cross-sectional volatility
avevol_2 = zeros(nsim, 1);
end

[em, sigm, Sp, lambda] = lambda_m(x, beta, gamA, gamB, xbar, rhox, stdx);
rf = 1./em;

for is   = 1 : nsim
    %--------------------------------------------------
    randn('state', sum(100*clock))

    % aggregate shocks 
    [sx, xinlnew] = CspSimu(xinl, xbar, rhox, stdx, Ts);
    xinl      = xinlnew;
    % stochastic discount factor
    SDF       = beta*exp((gamA + gamB*(sx(1:end-1) - xbar)).*(sx(1:end-1) - sx(2:end)));
    % real interest rate
    srf       = interp1(x, rf, sx, 'linear', 'extrap');    
    srf(end)  = [];
    % price of risk
    slam      = interp1(x, lambda, sx, 'linear', 'extrap');    
    slam(end) = [];
    
    % simulation MEX routine
    % [Pf, Bf, Df, Rf, In, iyr(is), ikr(is), theta(is), dyr(is), fyr(is), Rm, GDPg, kd0new, zd0new] = ...
    %    panIEfcn(kd0, zd0, optK, V0, k, sx, h, nh, x, nx, z, N, Ts, alpha, alp1, alp2, alp3, delta, eta, f, gP, gN, istar, rhoz, stdz);
    [Pf, Bf, Df, Rf, In, iyr(is), ikr(is), theta(is), dyr(is), fyr(is), Rm, GDPg, kd0new, zd0new, Zf] = ...
       panIEfcnZ(kd0, zd0, optK, V0, k, sx, h, nh, x, nx, z, N, Ts, alpha, alp1, alp2, alp3, delta, eta, f, gP, gN, istar, rhoz, stdz);
    
    Pf(:, end) = [];     Bf(:, end) = [];     Df(:, end) = [];    In(:, end) = [];    Zf(:, end) = [];
    
    % average rate of investment and disinvestment
    negInd     = find(In < 0);
    posInd     = find(In > 0);
    negI(is)   = mean(In(negInd) ./Bf(negInd));
    posI(is)   = mean(In(posInd) ./Bf(posInd));
    clear negInd posInd
    
    % initial distribution for next iteration
    kd0 = kd0new;
    zd0 = zd0new;
    % 
    % quality control
    if min(min(Pf)) < 0
        fprintf('negative firm value!\n');
        Pf = max(Pf, 1e-8);
    end
    
    sx(end) = [];
    
    % TIME SERIES ANALYSIS
    % 
    % [HML, mktrf] = bpDriver(Pf, Bf, Rf, Zf, Rm, srf, sx); 
    % sHML  = [sHML HML];
    % sMKT  = [sMKT Rm];
    % simrf = [simrf srf];

    % Pontiff-Schall
    btm   = Bf./Pf;
    sbtm  = sum(Bf)./sum(Pf); 
   
    meanbtm(:, is) = mean(mean(btm));
    volabtm(:, is) = mean(std(btm));

    btmPanel(:, :, is) = ponSchall(Rm(61 : end), sbtm(61 : end));
     
    % Industry return on the value spread
    [Junk, Junk0, ValSpr, Junk1, Junk2, Junk3] = p3by2SortEx(Pf, Rf, Bf, Df);          clear Junk Junk0 Junk1 Junk2 Junk3
    rovPanel(:, :, is)  = RmOnVp(Rm', ValSpr');
    % 
    % Equilibrium effects
    sigBM = std(btm);   
    sigR  = std(Rf);
    equiPanel(:, :, is) = EquiEffects(Rm', sx' - xbar, sigBM' - mean(sigBM), sigR' - mean(sigR));
    
    % CROSS-SECTIONAL ANALYSIS
    % 
    if 1
    % summary statistics of book value
    [mvec(:, is), bvec(:, is), sigvec(:, is)] = book10Sort(Pf, Bf, Rf, Zf, Rm, srf);
        
    [p10btm, k10btm, z10btm, p10size, k10size, z10size] = p10Sort(Pf, Bf, Rf, Zf);
    [HML, SMB] = p3by2Sort(Pf, Rf, Bf);
    % 
    % summary statistics of 10 portfolios
    [mbtm(:, is), stdbtm(:, is), msize(:, is), stdsize(:, is), mhml(is), vhml(is), msmb(is), vsmb(is), ...
            alp_size(:, is), talp_size(:, is), bet_size(:, is), tbet_size(:, is), ...
            alp_btm(:, is), talp_btm(:, is), bet_btm(:, is), tbet_btm(:, is), ...
            alp_hml(is), talp_hml(is), bet_hml(is), tbet_hml(is), ...
            alp_smb(is), talp_smb(is), bet_smb(is), tbet_smb(is), ...
            percbtm(:, is), percsize(:, is), perchml(is), percsmb(is)] = SummStats(p10btm, p10size, HML, SMB, Rm, srf);
    tmp = corrcoef(HML, SMB);
    corrhmlsmb(is) = tmp(1, 2);
    % 
    % [mbtm(:, is), stdbtm(:, is), msize(:, is), stdsize(:, is), ...
    %       mhml(is), vhml(is), msmb(is), vsmb(is), ...  % summary statistics
    %       alp_size(:, is), talp_size(:, is), bet_size(:, is), tbet_size(:, is), ...
    %       alp_btm(:, is), talp_btm(:, is), bet_btm(:, is), tbet_btm(:, is), ... % lapm regressions
    %       alp_hml(is), talp_hml(is), bet_hml(is), tbet_hml(is), alp_smb(is), talp_smb(is), bet_smb(is), tbet_smb(is), ...
    %       percbtm(:, is), percsize(:, is), perchml(is), percsmb(is)] = lapm(p10btm, p10size, HML, SMB, slam, srf);
    
    % beta-premium analysis
    % [out_btm(:, :, is)] = betaPrem(p10btm, k10btm, z10btm, HML, SMB, Rm, srf, sx);
    end
    % 
    % Fama and French 1995 Analysis --- need Df data from simulation
    % [roevec(:, :, is), btmvec(:, :, is), ROElow(is, :), ROEhigh(is, :)] = FF95(Pf, Bf, Df);

    if 0
    % value factor in investments
    % 
    % generate the firm index of growth and value at every point of time
    [GrowthId, ValueId] = p10SortIndex(Pf, Bf);
    % panels of adjustment cost and i/k ratio
    Af     = (gP/2)*(In./Bf).^2.*Bf;
    Id     = find(In < 0);
    Af(Id) = (gN/2)*(In(Id)./Bf(Id)).^2.*Bf(Id);   clear Id
    ikr    = In./Bf;   clear In 
    % good and bad times
    boom   = find(sx > mean(sx) + 1*stdx/sqrt(1 - rhox^2));
    bust   = find(sx < mean(sx) - 1*stdx/sqrt(1 - rhox^2));
    t_boom = boom(unidrnd(length(boom)));
    t_bust = bust(unidrnd(length(bust)));
    % 
    xgboom = ikr(GrowthId(:, t_boom), t_boom);
    ygboom = Af(GrowthId(:, t_boom), t_boom);
    xvboom = ikr(ValueId(:, t_boom), t_boom);
    yvboom = Af(ValueId(:, t_boom), t_boom);
    xgbust = ikr(GrowthId(:, t_bust), t_bust);
    ygbust = Af(GrowthId(:, t_bust), t_bust);
    xvbust = ikr(ValueId(:, t_bust), t_bust);
    yvbust = Af(ValueId(:, t_bust), t_bust);
    end

    % style timing
    [c_a(:, is), t_a(:, is), R2_a(:, is), corrmtrx_a(:, :, is)] = StyleTiming(Pf, Rf, Bf, Df, sx, 1);
    [c_m(:, is), t_m(:, is), R2_m(:, is), corrmtrx_m(:, :, is)] = StyleTiming(Pf, Rf, Bf, Df, sx, 0);

    fprintf('is = %6.0f   time = %6.4f\n', [is toc])
        
    % average volatility of individual stocks
    % avevol_1(is) = mean(std(Rf, 0, 2))*sqrt(12);
    % avevol_2(is) = mean(std(Rf))*sqrt(12);
    
%-----------------------------------------------------------
end  % nsim

posI = mean(posI)*12
negI = mean(negI)*12

if 1
mmbtm   = mean(mbtm, 2)
msbtm   = mean(stdbtm, 2)
mbbtm   = mean(bet_btm, 2)

mmhml   = mean(mhml)
mvhml   = mean(vhml)
mbhml   = mean(bet_hml)
end

meanbtm = mean(meanbtm)
volabtm = mean(volabtm)

% mvol_1  = mean(avevol_1)
% mvol_2  = mean(avevol_2)

c_m        = mean(c_m, 2)
t_m        = mean(t_m, 2)
R2_m       = mean(R2_m, 2)
corrmtrx_m = mean(corrmtrx_m, 3)

c_a        = mean(c_a, 2)
t_a        = mean(t_a, 2)
R2_a       = mean(R2_a, 2)
corrmtrx_a = mean(corrmtrx_a, 3)


btmmean    = mean(btmPanel, 3)
rovmean    = mean(rovPanel, 3)
equimean   = mean(equiPanel, 3)


return


% now plot!
figure(1) % boom
set(gca, 'FontS', 15);
plot(xgboom, ygboom, 'o'); hold on;
plot(xvboom, yvboom, 'r+'); hold off;
xlabel('i/k', 'FontS', 20); ylabel('Adjustment Cost', 'FontS', 20);
gtext('Growth', 'FontS', 20); gtext('Value', 'FontS', 20);
set(gca, 'FontS', 15);
print -deps scatter_boom.eps

figure(2) % bust
set(gca, 'FontS', 15);
plot(xgbust, ygbust, 'o'); hold on;
plot(xvbust, yvbust, 'r+'); hold off;
xlabel('i/k', 'FontS', 20); ylabel('Adjustment Cost', 'FontS', 20);
gtext('Growth', 'FontS', 20); gtext('Value', 'FontS', 20);
set(gca, 'FontS', 15);
print -deps scatter_bust.eps


% ff95 plots
roevec  = mean(roevec, 3);
btmvec  = mean(btmvec, 3);
ROElow  = mean(ROElow);
ROEhigh = mean(ROEhigh);

xseries = -5 : 1: 5;

% plot profitability: Figure 1 in FF 95
figure(3)
plot(xseries, roevec(:, 1), '-', 'LineWidth', 2); hold on; 
plot(xseries, roevec(:, 2), '--', 'LineWidth', 2); hold off;
xlabel('Formation Year', 'FontS', 15); ylabel('Profitability', 'FontS', 15);
set(gca, 'FontS', 15);
axis([-6 6 -0.2 0.4])
gtext('Growth', 'FontS', 20)
gtext('Value', 'FontS', 20)
print -deps roe11.eps

% plot book-to-market: Figure 2 in FF 95
% figure(2)
% plot(xseries, btmvec(:, 1), '-', 'LineWidth', 2); hold on; 
% plot(xseries, btmvec(:, 2), '--', 'LineWidth', 2); hold off;
% xlabel('Formation Year', 'FontS', 15); ylabel('k/v^e', 'FontS', 15);
% set(gca, 'FontS', 15);
% axis([-6 6 0.35 .85])
% gtext('Growth', 'FontS', 20)
% gtext('Value', 'FontS', 20)
% print -deps btm11.eps

% plot profitability: Figure 3 in FF 95
figure(4);
plot(ROElow(2:end), '-', 'LineWidth', 2); hold on; 
plot(ROEhigh(2:end), '--', 'LineWidth', 2); hold off; 
xlabel('Time Series', 'FontS', 15); ylabel('Profitability', 'FontS', 15);
set(gca, 'FontS', 15);
gtext('Growth', 'FontS', 20)
gtext('Value', 'FontS', 20)
print -deps profitability.eps

return



mmvec   = mean(mvec, 2)
mbvec   = mean(bvec, 2)
msigvec = mean(sigvec, 2)

mmsize   = mean(msize, 2)
mssize   = mean(stdsize, 2)
mbsize   = mean(bet_size, 2)

return

% report results
% hml
% mhml = mean(sHML) * 12
% stdh = std(sHML) * sqrt(12)
% 
% industry return
% mind = mean(Rm)^12
% vind = std(Rm)*sqrt(12)
% 
% real interest rate
meanrf = mean(simrf)^12
stdrf  = std(simrf)*sqrt(12)
% 
% industry btm
% btma   = TimeAggr(sbtm, 'm', 'y', 's');
% mbtm   = mean(btma)
% vbtm   = std(btma)
% 
% ratios
% mtheta = mean(theta)
% miyr   = mean(iyr)
% mpikr  = mean(posI)*12
% mnikr  = mean(negI)*12
% 
% btmPanMean = mean(btmPanel, 3)
% btmPanStd  = std(btmPanel, 0, 3)

if 0
mc    = mean(c, 2)
mt    = mean(t, 2)
mR2   = mean(R2, 2)
mcorr = mean(corrmtrx, 3)
end


if 0
mmbtm   = mean(mbtm, 2)
msbtm   = mean(stdbtm, 2)
mmsize  = mean(msize, 2)
mssize  = mean(stdsize, 2)
mmhml   = mean(mhml)
mvhml   = mean(vhml)
mmsmb   = mean(msmb)
mvsmb   = mean(vsmb)
masize  = mean(alp_size, 2)
mtasize = mean(talp_size, 2)
mbsize  = mean(bet_size, 2)
mtbsize = mean(tbet_size, 2)
mabtm   = mean(alp_btm, 2)
mtabtm  = mean(talp_btm, 2)
mbbtm   = mean(bet_btm, 2)
mtbbtm  = mean(tbet_btm, 2)
mahml   = mean(alp_hml)
mtahml  = mean(talp_hml)
mbhml   = mean(bet_hml)
mtbhml  = mean(tbet_hml)
masmb   = mean(alp_smb)
mtasmb  = mean(talp_smb)
mbsmb   = mean(bet_smb)
mtbsmb  = mean(tbet_smb)
mmpercbtm  = mean(percbtm, 2)
mmpercsize = mean(percsize, 2)
mmperchml  = mean(perchml)
mmpercsmb  = mean(percsmb)
end

% 
% rovPan1Mean = mean(rovPanel1, 3)
% rovPan2Mean = mean(rovPanel2, 3)
% 
% save artipanel Pf Bf Rf Df In Zf Rm srf sx

return

Xm = [ones(size(sMKT')) sMKT'];
[c, t, R2, chi2, p] = MultiRegressNw(sHML', Xm, [], [], 0);
alp_hml = c(1); talp_hml = t(1);
bet_hml = c(2); tbet_hml = t(2);

% some simple summary statistics
% 
% equity premium
meanEP = mean(simEP)*12
stdEP  = std(simEP)*sqrt(12)
% 
% 
% btm
mBTM   = mean(simBTM)
sBTM   = std(simBTM)
% 
% aggregate ratios
shares = mean([iyr ikr theta dyr fyr negI posI])