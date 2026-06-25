
function [mbtm, stdbtm, msize, stdsize, mhml, vhml, msmb, vsmb, ...
          alp_size, talp_size, bet_size, tbet_size, alp_btm, talp_btm, bet_btm, tbet_btm, ...
          alp_hml, talp_hml, bet_hml, tbet_hml, alp_smb, talp_smb, bet_smb, tbet_smb, ...
          percbtm, percsize, perchml, percsmb] = SummStats(p10btm, p10size, hml, smb, Rm, srf)

% -------------------------------------------------------------------
% This file reports the summary statistics and unconditional market
% regression in the beta-premium analysis
% -------------------------------------------------------------------

% mean and volatility
mbtm    = mean(p10btm, 2).^12;
stdbtm  = std(p10btm, 0, 2)*sqrt(12);
msize   = mean(p10size, 2).^12;
stdsize = std(p10size, 0, 2)*sqrt(12);
mhml    = mean(hml)*12;
vhml    = std(hml)*sqrt(12);
msmb    = mean(smb)*12;
vsmb    = std(smb)*sqrt(12);

mktrf   = (Rm - srf)';    
p10btm  = (p10btm - repmat(srf, 10, 1))';
p10size = (p10size - repmat(srf, 10, 1))';
%
Xm = [ones(size(mktrf)) mktrf];
alp_size = zeros(10, 1);  talp_size = zeros(10, 1); 
bet_size = zeros(10, 1);  tbet_size = zeros(10, 1);
alp_btm  = zeros(10, 1);  talp_btm  = zeros(10, 1);
bet_btm  = zeros(10, 1);  tbet_btm  = zeros(10, 1);
% 
for j = 1 : 10
   % btm portfolio
    port = p10btm(:, j);
    [c, t, R2, chi2, p] = MultiRegressNw(port, Xm, [], [], 0);
    alp_btm(j) = c(1); talp_btm(j) = t(1);
    bet_btm(j) = c(2); tbet_btm(j) = t(2);
    % size portfolio
    port = p10size(:, j);
    [c, t, R2, chi2, p] = MultiRegressNw(port, Xm, [], [], 0);
    alp_size(j) = c(1); talp_size(j) = t(1);
    bet_size(j) = c(2); tbet_size(j) = t(2);
end

mmktrf = mean(mktrf);
% 
percbtm  = (bet_btm*mmktrf) ./mean(p10btm)';
percsize = (bet_size*mmktrf)./mean(p10size)';

% hml
[c, t, R2, chi2, p] = MultiRegressNw(hml', Xm, [], [], 0);
alp_hml = c(1); talp_hml = t(1);
bet_hml = c(2); tbet_hml = t(2);
% smb
[c, t, R2, chi2, p] = MultiRegressNw(smb', Xm, [], [], 0);
alp_smb = c(1); talp_smb = t(1);
bet_smb = c(2); tbet_smb = t(2);

%
perchml = bet_hml*mmktrf/mean(hml);
percsmb = bet_smb*mmktrf/mean(smb);