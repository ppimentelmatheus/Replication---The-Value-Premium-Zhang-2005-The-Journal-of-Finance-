
function [roevec, btmvec, ROElow, ROEhigh] = FF95(V_f, B_f, D_f)

% FF95:
%   Examine the behavior of profitability measured by return on book equity (ROEs):
% 
%                    ROE_{t+1} = log((BE_{t+1} + D_t)/BE_t)
% 
% and its relation to BE/ME 
% 
% OUTPUT:
%   roevec  -- 11 by 2 matrix
%   btmvec  -- 11 by 2 matrix
%   ROElow  -- 1 by T/12 (= 35 in my simulation) matrix
%   ROEhigh -- 1 by T/12 matrix
% 
% NOTES: Fama and French (1995)
%
% © Lu Zhang, Inc. 2001


% initializing
[N, T] = size(V_f);
% profitability of FF 95
ROE    = (B_f + D_f) ./east(B_f) - 1;
% earnings
E      = B_f + D_f - east(B_f);
btm    = B_f ./V_f;                

% sindex range of top and bottom decile returns
gpSize = round(.30*N);
low    = 1 : gpSize;
high   = N - gpSize + 1 : N;

% 11-year evoluation of profitability and B/M: FF 95 Figure 1 and 2
roeM5 = [];  roeM4 = [];  roeM3 = [];  roeM2 = [];  roeM1 = [];  
roe00 = [];
roeP1 = [];  roeP2 = [];  roeP3 = [];  roeP4 = [];  roeP5 = [];
% 
btmM5 = [];  btmM4 = [];  btmM3 = [];  btmM2 = [];  btmM1 = [];
btm00 = [];
btmP1 = [];  btmP2 = [];  btmP3 = [];  btmP4 = [];  btmP5 = [];

for t = 62 : T - 71
    if mod(t - 2, 12) == 0,
       [Junk, sindex]  = sort(btm(:, t));        clear Junk;
       % 11-year evoluation of profitability
       roeM5 = [roeM5; 
           12*mean(sum(E(sindex(low),  t - 60 : t - 49))./sum(B_f(sindex(low),  t - 61 : t - 50))) ...
           12*mean(sum(E(sindex(high), t - 60 : t - 49))./sum(B_f(sindex(high), t - 61 : t - 50))) ];
       roeM4 = [roeM4; 
           12*mean(sum(E(sindex(low),  t - 48 : t - 37))./sum(B_f(sindex(low),  t - 49 : t - 38))) ...
           12*mean(sum(E(sindex(high), t - 48 : t - 37))./sum(B_f(sindex(high), t - 49 : t - 38))) ];
       roeM3 = [roeM3; 
           12*mean(sum(E(sindex(low),  t - 36 : t - 25))./sum(B_f(sindex(low),  t - 37 : t - 26))) ...
           12*mean(sum(E(sindex(high), t - 36 : t - 25))./sum(B_f(sindex(high), t - 37 : t - 26))) ];
       roeM2 = [roeM2; 
           12*mean(sum(E(sindex(low),  t - 24 : t - 13))./sum(B_f(sindex(low),  t - 25 : t - 14))) ...
           12*mean(sum(E(sindex(high), t - 24 : t - 13))./sum(B_f(sindex(high), t - 25 : t - 14))) ];
       roeM1 = [roeM1; 
           12*mean(sum(E(sindex(low),  t - 12 : t - 1))./sum(B_f(sindex(low),  t - 13 : t - 2))) ...
           12*mean(sum(E(sindex(high), t - 12 : t - 1))./sum(B_f(sindex(high), t - 13 : t - 2))) ];
       roe00 = [roe00; 
           12*mean(sum(E(sindex(low),  t : t + 11))./sum(B_f(sindex(low),  t - 1 : t + 10))) ...
           12*mean(sum(E(sindex(high), t : t + 11))./sum(B_f(sindex(high), t - 1 : t + 10))) ];
       roeP1 = [roeP1; 
           12*mean(sum(E(sindex(low),  t + 12 : t + 23))./sum(B_f(sindex(low),  t + 11 : t + 22))) ...
           12*mean(sum(E(sindex(high), t + 12 : t + 23))./sum(B_f(sindex(high), t + 11 : t + 22))) ];
       roeP2 = [roeP2; 
           12*mean(sum(E(sindex(low),  t + 24 : t + 35))./sum(B_f(sindex(low),  t + 23 : t + 34))) ...
           12*mean(sum(E(sindex(high), t + 24 : t + 35))./sum(B_f(sindex(high), t + 23 : t + 34))) ];
       roeP3 = [roeP3; 
           12*mean(sum(E(sindex(low),  t + 36 : t + 47))./sum(B_f(sindex(low),  t + 35 : t + 46))) ...
           12*mean(sum(E(sindex(high), t + 36 : t + 47))./sum(B_f(sindex(high), t + 35 : t + 46))) ];
       roeP4 = [roeP4; 
           12*mean(sum(E(sindex(low),  t + 48 : t + 59))./sum(B_f(sindex(low),  t + 47 : t + 58))) ...
           12*mean(sum(E(sindex(high), t + 48 : t + 59))./sum(B_f(sindex(high), t + 47 : t + 58))) ];
       roeP5 = [roeP5; 
           12*mean(sum(E(sindex(low),  t + 60 : t + 71))./sum(B_f(sindex(low),  t + 59 : t + 70))) ...
           12*mean(sum(E(sindex(high), t + 60 : t + 71))./sum(B_f(sindex(high), t + 59 : t + 70))) ];
    end
end

roevec = [  mean(roeM5);
            mean(roeM4);
            mean(roeM3);
            mean(roeM2);
            mean(roeM1);
            mean(roe00);
            mean(roeP1);
            mean(roeP2);
            mean(roeP3);
            mean(roeP4);
            mean(roeP5) ];

for t = 61 : T - 71
    if mod(t - 1, 12) == 0,
       [Junk, sindex]  = sort(btm(:, t));        clear Junk;
       % 11-year evoluation of btm
       btmM5 = [btmM5; mean(mean(btm(sindex(low), t - 60 : t - 49)))  mean(mean(btm(sindex(high), t - 60 : t - 49)))];
       btmM4 = [btmM4; mean(mean(btm(sindex(low), t - 48 : t - 37)))  mean(mean(btm(sindex(high), t - 48 : t - 37)))];
       btmM3 = [btmM3; mean(mean(btm(sindex(low), t - 36 : t - 25)))  mean(mean(btm(sindex(high), t - 36 : t - 25)))];
       btmM2 = [btmM2; mean(mean(btm(sindex(low), t - 24 : t - 13)))  mean(mean(btm(sindex(high), t - 24 : t - 13)))];
       btmM1 = [btmM1; mean(mean(btm(sindex(low), t - 12 : t - 1)))   mean(mean(btm(sindex(high),  t - 12 : t - 1)))];
       btm00 = [btm00; mean(mean(btm(sindex(low), t : t + 11)))       mean(mean(btm(sindex(high), t : t + 11)))     ];
       btmP1 = [btmP1; mean(mean(btm(sindex(low), t + 12 : t + 23)))  mean(mean(btm(sindex(high), t + 12 : t + 23)))];
       btmP2 = [btmP2; mean(mean(btm(sindex(low), t + 24 : t + 35)))  mean(mean(btm(sindex(high), t + 24 : t + 35)))];
       btmP3 = [btmP3; mean(mean(btm(sindex(low), t + 36 : t + 47)))  mean(mean(btm(sindex(high), t + 36 : t + 47)))];
       btmP4 = [btmP4; mean(mean(btm(sindex(low), t + 48 : t + 59)))  mean(mean(btm(sindex(high), t + 48 : t + 59)))];
       btmP5 = [btmP5; mean(mean(btm(sindex(low), t + 60 : t + 71)))  mean(mean(btm(sindex(high), t + 60 : t + 71)))];
   end
end
        
btmvec = [  mean(btmM5);
            mean(btmM4);
            mean(btmM3);
            mean(btmM2);
            mean(btmM1);
            mean(btm00);
            mean(btmP1);
            mean(btmP2);
            mean(btmP3);
            mean(btmP4);
            mean(btmP5) ];

% Profitability
ROElow    = zeros(1, T);
ROEhigh   = zeros(1, T);
% Earnings
% Elow      = zeros(1, T);
% Ehigh     = zeros(1, T);

% Time Series Evidence: FF 95 Figure 3
for t = 1 : T - 11
   % Construct value premium based on sorting on btm
   if mod(t-1, 12) ==  0,
       [Junk, sindex]  = sort(btm(:, t));        clear Junk;
       % Profitability
       ROElow(t : t + 11)  = mean(ROE(sindex(low), t : t + 11));
       ROEhigh(t : t + 11) = mean(ROE(sindex(high), t : t + 11));
       % Earnings
       % Elow(t : t + 11)    = mean(E(sindex(low), t : t + 11));
       % Ehigh(t : t + 11)   = mean(E(sindex(high), t : t + 11));
   end
end; clear sindex

% Time Aggregation to Annual Frequency
if 1, % FF 95
    ROElow  = 12*mean(reshape(ROElow, 12, T/12));
    ROEhigh = 12*mean(reshape(ROEhigh, 12, T/12));
else, % Cohen, Polk, and Vuolteenaho 2000
    ROElow  = sum(reshape(ROElow, 12, T/12));
    ROEhigh = sum(reshape(ROEhigh, 12, T/12));
end    










return

% Some coding junk follows ---
% 
% Plot earnings: levels
Elow   = sum(reshape(Elow, 12, T/12));
Ehigh  = sum(reshape(Ehigh, 12, T/12));
figure(2);
plot(Elow, '-', 'LineWidth', 2); hold on;
plot(Ehigh, ':', 'LineWidth', 2); hold off; 
xlabel('Time Series', 'FontS', 15); ylabel('Earnings Level', 'FontS', 15);
set(gca, 'FontS', 20);
legend('Low B/M', 'High B/M');
print -deps e:\Research\ValPrem\Documents\FIGURES\earnLevel.eps

% Plot earnings: growth rate
EGlow  = Elow(2 : end)./Elow(1 : end - 1);
EGhigh = Ehigh(2 : end)./Ehigh(1 : end - 1);
figure(3);
plot(EGlow(2:end), '-', 'LineWidth', 2); hold on;
plot(EGhigh(2:end), ':', 'LineWidth', 2); hold off; 
xlabel('Time Series', 'FontS', 15); ylabel('Earnings Growth', 'FontS', 15);
set(gca, 'FontS', 20);
legend('Low B/M', 'High B/M');
print -deps e:\Research\ValPrem\Documents\FIGURES\earnGrowth.eps

% Plot time series of ROEs for two typical firms
roe1firm   = 12*mean(reshape(ROE(100,  :), 12, T/12));
roe2firm   = 12*mean(reshape(ROE(1000, :), 12, T/12));
figure(4);
plot(roe1firm, '-', 'LineWidth', 2); hold on;
plot(roe2firm, ':', 'LineWidth', 2); hold off; 
xlabel('Time Series', 'FontS', 15); ylabel('ROE', 'FontS', 15);
set(gca, 'FontS', 20);
legend('Firm 1', 'Firm 2');
print -deps e:\Research\ValPrem\Documents\FIGURES\roeFirms.eps

% Average First-order autocorrelation
autocorr = zeros(N_f, 1);
for f = 1 : N_f
    roefirm = 12*mean(reshape(ROE(f, :), 12, T/12));
    tmp     = corrcoef(roefirm(1:end-1), roefirm(2:end));
    autocorr(f) = tmp(1, 2);
end

mean(autocorr)

return

% some coding junk
if 0 % A different approach for computing accounting ratios for portfolios
for t = 61 : T - 71
    if mod(t - 1, 12) == 0,
       [Junk, sindex]  = sort(btm(:, t));        clear Junk;
       % 11-year evoluation of profitability
       roeM5 = [roeM5; 12*mean(mean(ROE(sindex(low), t - 60 : t - 49))) 12*mean(mean(ROE(sindex(high), t - 60 : t - 49)))];
       roeM4 = [roeM4; 12*mean(mean(ROE(sindex(low), t - 48 : t - 37))) 12*mean(mean(ROE(sindex(high), t - 48 : t - 37)))];
       roeM3 = [roeM3; 12*mean(mean(ROE(sindex(low), t - 36 : t - 25))) 12*mean(mean(ROE(sindex(high), t - 36 : t - 25)))];
       roeM2 = [roeM2; 12*mean(mean(ROE(sindex(low), t - 24 : t - 13))) 12*mean(mean(ROE(sindex(high), t - 24 : t - 13)))];
       roeM1 = [roeM1; 12*mean(mean(ROE(sindex(low), t - 12 : t - 1)))  12*mean(mean(ROE(sindex(high),  t - 12 : t - 1)))];
       roe00 = [roe00; 12*mean(mean(ROE(sindex(low), t : t + 11)))      12*mean(mean(ROE(sindex(high), t : t + 11)))];
       roeP1 = [roeP1; 12*mean(mean(ROE(sindex(low), t + 12 : t + 23))) 12*mean(mean(ROE(sindex(high), t + 12 : t + 23)))];
       roeP2 = [roeP2; 12*mean(mean(ROE(sindex(low), t + 24 : t + 35))) 12*mean(mean(ROE(sindex(high), t + 24 : t + 35)))];
       roeP3 = [roeP3; 12*mean(mean(ROE(sindex(low), t + 36 : t + 47))) 12*mean(mean(ROE(sindex(high), t + 36 : t + 47)))];
       roeP4 = [roeP4; 12*mean(mean(ROE(sindex(low), t + 48 : t + 59))) 12*mean(mean(ROE(sindex(high), t + 48 : t + 59)))];
       roeP5 = [roeP5; 12*mean(mean(ROE(sindex(low), t + 60 : t + 71))) 12*mean(mean(ROE(sindex(high), t + 60 : t + 71)))];
       % 11-year evoluation of profitability
       btmM5 = [btmM5; mean(mean(btm(sindex(low), t - 60 : t - 49)))  mean(mean(btm(sindex(high), t - 60 : t - 49)))];
       btmM4 = [btmM4; mean(mean(btm(sindex(low), t - 48 : t - 37)))  mean(mean(btm(sindex(high), t - 48 : t - 37)))];
       btmM3 = [btmM3; mean(mean(btm(sindex(low), t - 36 : t - 25)))  mean(mean(btm(sindex(high), t - 36 : t - 25)))];
       btmM2 = [btmM2; mean(mean(btm(sindex(low), t - 24 : t - 13)))  mean(mean(btm(sindex(high), t - 24 : t - 13)))];
       btmM1 = [btmM1; mean(mean(btm(sindex(low), t - 12 : t - 1)))   mean(mean(btm(sindex(high),  t - 12 : t - 1)))];
       btm00 = [btm00; mean(mean(btm(sindex(low), t : t + 11)))       mean(mean(btm(sindex(high), t : t + 11)))     ];
       btmP1 = [btmP1; mean(mean(btm(sindex(low), t + 12 : t + 23)))  mean(mean(btm(sindex(high), t + 12 : t + 23)))];
       btmP2 = [btmP2; mean(mean(btm(sindex(low), t + 24 : t + 35)))  mean(mean(btm(sindex(high), t + 24 : t + 35)))];
       btmP3 = [btmP3; mean(mean(btm(sindex(low), t + 36 : t + 47)))  mean(mean(btm(sindex(high), t + 36 : t + 47)))];
       btmP4 = [btmP4; mean(mean(btm(sindex(low), t + 48 : t + 59)))  mean(mean(btm(sindex(high), t + 48 : t + 59)))];
       btmP5 = [btmP5; mean(mean(btm(sindex(low), t + 60 : t + 71)))  mean(mean(btm(sindex(high), t + 60 : t + 71)))];
   end
end
end
