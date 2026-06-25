% MAINCC:
%   Construct the industry equilibrium without entry and exit by applying Krusell and Smith [1998] method
% 
%                  p' = alp1 + alp2 p + alp3 (x - xbar)
%      log(M_{t,t+1}) = log(beta) + (gamA + gamB*(x_t - xbar))*(x_t - x_{t+1})
% 
% MODIFICATION LOG:
%   08/22/01 -- rewrite IE to (1) call vfi3fcn and simulation procedures as functional files and (2) use 
%        continuous state space method in the simulation (3) vfi3fcn and simulation will be MEX files
%   12/22/01 -- implement Campbell and Cochrane pricing kernel
% 
% © Lu Zhang, Inc. 2001

clear all; clc; format compact; format short; 

load Params

if 1
    load coefIS
else
    % initializing approximating coefficients --- log average capital stock is linear: h' = alp1 + alp2 h + alp3 (x - xbar)
    alp1  = 0;         
    alp2  = 1;         
    alp3  = 0;      
    alp4  = 0;
end
% initializing the economy
cutoff = 1000;
xinl   = xbar;   
% 
% aggregate shocks on discrete state space
sx    = CspSimu(xinl, xbar, rhox, stdx, Ts);
% simulated path of aggregate shock to be used in regressions
simx  = sx(cutoff + 1 : end);
% 
% initializing 
kd    = ones(N, 1);                                       % cross-sectional current period capital stock
% idiosyncratic shocks
shock = randn(N/2, 1);                                    % impose law of large numbers using antithetic methods
zd    = (stdz/sqrt(1 - rhoz^2)) * [shock; -shock];        % idiosyncratic shock from continuous state space

% some intermediate variables used in vfi3fcn routine -- follow Qzxh's convention of rate of change in the second dimension
kmtrx = repmat(k, 1, nz*nx*nh);
ksubm = repmat(k, 1, nz*nh);
hmtrx = repmat(repmat(h, nz*nx, 1)', nk, 1);                              % h changes the fastest
xmtrx = repmat(repmat(kron(x, ones(nh, 1)), nz, 1)', nk, 1);              % x changes the second fastest
zmtrx = repmat(kron(z, ones(nx*nh, 1))', nk, 1);                          % z changes the slowest

if 1, load vfi3Mat; else, V0 = kmtrx; end

% starting krusell and smith iteration 
err  = 1;
iter = 0;

tic;
while err >= 1e-2
    iter = iter + 1;
  
    % for given approximating parameters solve for the value function and optimal decision rules
    [optK, V, I, div] = vfi3fcnIEccB(alp1, alp2, alp3, V0, k, nk, x, xbar, nx, Qx, z, nz, Qz, h, nh, kp, ...
                                  alpha, beta, delta, f, gamA, gamB, gP, gN, istar, kmin, kmtrx, ksubm, hmtrx, xmtrx, zmtrx);
    % for given approximating parameters solve for the value function and optimal decision rules
    % [optK, V, Va, I, div] = vafcnIEcc(alp1, alp2, alp3, V0, k, nk, x, xbar, nx, Qx, z, nz, Qz, h, nh, kp, ...
    %                             alpha, beta, delta, f, gamA, gamB, gP, gN, istar, kmin, kmtrx, ksubm, hmtrx, xmtrx, zmtrx);
    V0 = V;         
    save vfi3Mat V0 optK
    save PD V optK I div k x z h 
    
    % reshape to facilitate interpolation later
    optKtmp = reshape(optK, [nk nh nx nz]);
    
    % update the approximating coefficients by simulating the economy over a long time period
    % [simh] = simIEfcn(kd, zd, sx, optKtmp, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz);
    % [simh, sigk] = simIEfcn2(kd, zd, sx, optKtmp, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz);
    [simh, sigk, sigz] = simIEfcn3(kd, zd, sx, optKtmp, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz);
 
    % discard first 1000 periods to avoid inference of initial conditions of kd
    simh(1 : cutoff) = [];
    sigk(1 : cutoff) = [];    
    sigz(1 : cutoff) = [];    
    
    % update approximating coefficients on log average capital stock: h' = alp1 + alp2 h + alp3 (x - xbar) 
    % [coef, tstats, R2p, sighat, resid] = ols(simh(2:end)', [ones(length(simh) - 1, 1) simh(1:end-1)' simx(1:end-1)' - xbar]);
    % [coef, tstats, R2p, sighat, resid] = ols(simh(2:end)', [ones(length(simh) - 1, 1) simh(1:end-1)' ...
    %      simx(1 : end-1)' - xbar sigk(1 : end-1)' sigz(1 : end-1)']);
    [coef, tstats, R2p, sighat, resid] = ols(simh(2:end)', [ones(length(simh) - 1, 1) simh(1:end-1)' ...
           simx(1 : end-1)' - xbar sigk(1 : end-1)']);
    alp1new = coef(1);     
    alp2new = coef(2);     
    alp3new = coef(3);
    alp4new = coef(4);     
    % alp5new = coef(5);
    % convergence check
    % err   = max(abs([alp1new alp2new alp3new] - [alp1 alp2 alp3]));
    err  = max(abs([alp1new alp2new alp3new alp4new] - [alp1 alp2 alp3 alp4]));
    % prepare for next iteration
    alp1 = alp1new;    
    alp2 = alp2new;    
    alp3 = alp3new;
    alp4 = alp4new;    
    % alp5 = alp5new;
    
    % Visual
    fprintf('err = %6.4f    elapsed time: %6.4f\n', [err toc]);    
    save coefIS alp1 alp2 alp3 alp4 R2p sighat tstats resid simh
end

fprintf('competitive equilibrium is successfully constructed --- good job!\n')

% report results
coefs = [alp1 alp2 alp3 alp4]
tstats
R2p
sighat

% hist(resid'./simh(2 : end) *100, 100)
% xlabel('Percentage Deviation of Predicted Price from Actual Price', 'FontS', 20); ylabel('Frequency', 'FontS', 20);
% set(gca, 'FontS', 15);

predh = alp1 + alp2*simh(1 : end - 1) + alp3*(simx(1 : end - 1) - xbar) + alp4*sigk(1 : end - 1);

figure(1);
plot(predh/mean(predh), simh(2 : end)/mean(simh(2 : end)), '.', 'LineWidth', 2)
xlabel('Predicted Output Price', 'FontS', 20); ylabel('Actual Output Price', 'FontS', 20);
set(gca, 'FontS', 15); axis([.95 1.08 .95 1.08])
print -deps ..\Documents\FIGURES\predsimh.eps

figure(2);
hist(resid./simh(2 : end)' *100, 100)
xlabel('Excess Demand As a Precentage of Actual Output', 'FontS', 20); ylabel('Frequency', 'FontS', 20);
set(gca, 'FontS', 15);
print -deps ..\Documents\FIGURES\percdev.eps

ssIE