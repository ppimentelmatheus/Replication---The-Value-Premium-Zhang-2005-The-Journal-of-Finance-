
function [table1, table2, table3, table4] = shortPort(R_f, V_f, B_f, R, Portmtrx, postbeta, FirmIndex, N_f)

% SHORTPORT:
%   calculate summary statistics across portfolios
%
% OUTPUTS: 
%   table1 - 11 by 11 matrix of average returns for portfolios formed on size and pre-beta
%   table2 -  4 by 12 matrix of properties of portfolios formed on size
%   table3 -  4 by 12 matrix of properties of portfolios formed on prebeta
%   table4 -  4 by 12 matrix of properties of portfolios formed on book-to-market
%
% NOTES: table numbering follows that in fama and french [1992]
%   
% See also PORTFOLIOS for a more complete version from Gomes, Kogan, and Zhang (2001)
% 
% © Lu Zhang, Inc. 2000

% initialization
lag    = 60;
T      = size(R_f, 2);

% SIMPANEL simulates 3500+ firms as in ff [93] but ff [92] only use N_f firms
% randomly draw N_f firms as a sample from the panel with 3500+ firms
% vecsub = unidrnd(N, N_f, 1);
vecsub = 1 : 1 : N_f;
R_f    = R_f(vecsub, :);               
v_f    = log(V_f(vecsub, :));            clear V_f
btm    = log(B_f(vecsub, :)) - v_f;      clear B_f vecsub
% shift the timing of size and true beta to match the convention used in simulation
btm  = [zeros(N_f, 6)    btm(:, 1 : T - 6)];

% cut-off locations from ff (1992)
sizecutoff = [772+189 236 170 144 140 128 125 119 114 124   ];
betacutoff = [116 80 185 181 179 182 185 205 227 267 165 291];
btmcutoff  = [89 98 209 222 226 230 235 237 239 239 120 117 ];
% sizecutoff   = 200*ones(1, 10);
% betacutoff   = [100 100 200*ones(1, 8) 100 100];
% btmcutoff    = [100 100 200*ones(1, 8) 100 100];
numberfirms  = sum(sizecutoff);
numfirmsbeta = sum(betacutoff);
numfirmsbtm  = sum(btmcutoff);
for k = 2 : 10
   sizecutoff(k) = sizecutoff(k - 1) + sizecutoff(k); 
end
for k =2 : 12
   betacutoff(k) = betacutoff(k - 1) + betacutoff(k); 
   btmcutoff (k) = btmcutoff (k - 1) + btmcutoff (k);
end
sizecutoff = [1 1 + floor(sizecutoff./numberfirms *N_f)];
betacutoff = [1 1 + floor(betacutoff./numfirmsbeta*N_f)];
btmcutoff  = [1 1 + floor(btmcutoff ./numfirmsbtm *N_f)];
% cutoffs for 12 pure size portfolios
sizecutoff2 = [1 1 + floor((sizecutoff(2)-1)/2)  sizecutoff(2:10)...
      sizecutoff(10) + floor((sizecutoff(11)-sizecutoff(10))/2) sizecutoff(11) ];
Msize = zeros(12, N_f);
for j = 0 : 11
   Imin = sizecutoff2(1 + j);
   Imax = sizecutoff2(2 + j) - 1;
   Msize(j+1, Imin : Imax) = ones(1,Imax-Imin+1)./(Imax-Imin+1);  
end
Mbeta = zeros(12, N_f);
Mbtm  = zeros(12, N_f); 
for j = 0:11
   Imin = betacutoff(1+j);
   Imax = betacutoff(2+j)-1;
   Mbeta(j+1, Imin : Imax)  = ones(1, Imax-Imin+1)./(Imax-Imin+1);
   Imin = btmcutoff(1+j);
   Imax = btmcutoff(2+j)-1;
   Mbtm (j+1, Imin : Imax)  = ones(1, Imax-Imin+1)./(Imax-Imin+1);
end

% initialization
statsizeret  = zeros(12, T - lag - 1);
statsizebeta = zeros(12, T - lag - 1);
statsizesize = zeros(12, T - lag - 1);
statsizebtm  = zeros(12, T - lag - 1);
statbetaret  = zeros(12, T - lag - 1);
statbetabeta = zeros(12, T - lag - 1);
statbetasize = zeros(12, T - lag - 1);
statbetabtm  = zeros(12, T - lag - 1);
statbtmret   = zeros(12, T - lag - 1);
statbtmbeta  = zeros(12, T - lag - 1);
statbtmsize  = zeros(12, T - lag - 1);
statbtmbtm   = zeros(12, T - lag - 1);

for t = lag + 2 : T
   if (mod(t-1, 12) == 1)     
      % pre-beta estimation
      X       = [ones(1, lag); R(t - lag : t - 1); R(t - lag - 1 : t - 2)]';   
	  regcoef = inv(X'*X)*X'* R_f(:, t - lag : t - 1)';
	  prebeta = (regcoef(2, :) + regcoef(3, :))';
      % size and btm are recorded at the begining of the month
      sizevec = v_f(:, t);
      btmvec  = btm(:, t);
      % first sort by size
      Data      = [sizevec prebeta [1 : 1 : N_f]'];
      [y, i]    = sort(sizevec);	  clear y;
      databeta  = Data(:, 2);			  
      [z, n]    = sort(databeta);     clear z;	
      % sort by btm
      [z, nbtm] = sort(btmvec);       clear z;
      nsize     = i;
   end
   % size portfolios statistics
   statsizeret (:, t-lag-1) = Msize*R_f(nsize, t);
   Xtmp                     = [FirmIndex(:, t - lag - 1) postbeta];
   Xtmp                     = sortrows(Xtmp, 1);
   statsizebeta(:, t-lag-1) = Msize*Xtmp(nsize, 2);
   statsizesize(:, t-lag-1) = Msize*v_f(nsize,  t);
   statsizebtm (:, t-lag-1) = Msize*btm(nsize,  t);
   % beta portfolios statistics
   statbetaret (:, t-lag-1) = Mbeta*R_f (n, t);
   statbetabeta(:, t-lag-1) = Mbeta*Xtmp(n, 2);
   statbetasize(:, t-lag-1) = Mbeta*Data(n, 1);
   statbetabtm (:, t-lag-1) = Mbeta*btm (n, t);    
   % btm portfolios statistics
   statbtmret (:, t-lag-1)  = Mbtm*R_f (nbtm, t);
   statbtmbeta(:, t-lag-1)  = Mbtm*Xtmp(nbtm, 2);
   statbtmsize(:, t-lag-1)  = Mbtm*Data(nbtm, 1);
   statbtmbtm (:, t-lag-1)  = Mbtm*btm (nbtm, t); 
end

% average returns for portfolios formed on size and pre-beta
table1 = reshape(mean(Portmtrx, 2), 10, 10)';
table1 = [ mean(mean(table1)) mean(table1);
           mean(table1, 2)    table1       ] * 100; 

% post-ranking betas for portfolios formed on size and then pre-beta
% X              = [ones(1, T-lag-1); R(lag+1 : T-1); R(lag+2 : T)]';   
% regcoef        = inv(X'*X)*X'*Portmtrx';
% postbetaport   = (regcoef(2, :) + regcoef(3, :))';
% table5         = reshape(postbetaport, 10, 10)';
% table5         = [        0         mean(table5);
%                    mean(table5, 2)  table5        ];

% properties of portfolios formed on size
table2 = [ 100*mean(statsizeret,  2)';  mean(statsizebeta, 2)'; ...
               mean(statsizesize, 2)';  mean(statsizebtm,  2)' ];

% properties of portfolios formed on pre-beta
table3 = [ 100*mean(statbetaret,2)';    mean(statbetabeta,2)'; ...
               mean(statbetasize,2)';   mean(statbetabtm,2)' ];

% properties of portfolios formed on btm
table4 = [ 100*mean(statbtmret,  2)'; mean(statbtmbeta, 2)'; ...
               mean(statbtmsize, 2)'; mean(statbtmbtm,  2)' ];
