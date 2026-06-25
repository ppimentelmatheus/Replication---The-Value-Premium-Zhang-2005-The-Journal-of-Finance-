
function [ff, corr, Portmtrx, postbeta, FirmIndex] = FF92(V_f, R_f, B_f, R, N_f) 

% FF92:
%   perform on simulated panel three cross-sectional calculations:
%     (1) Fama and French [1992] cross-sectional regression
%     (2) cross-sectional correlations
%
% USAGE:
%   V_f    - panel of firm size
%   R_f    - panel of firm return
%   B_f    - panel of firm book value
%   R      - time series of market return
% 
% OUTPUT:
%   on cross-sectional analysis:
%   
%   ff     -  2 by 9 matrix containing results of fama-french regression:
%                   the first row is coefficients and the second row is t-statistics
%   corr   -  a 1 by 3 row vector of cross-sectional correlations
%   
% on intermediary matrix to be used in PORTFOLIOS:
% 
%      Portmtrx  -  matrix of 100 size-prebeta portfolios return series (100 by T)
%      postbeta  -  post-ranking beta of perfect ordering
%      FirmIndex -  size-prebeta double sorting index of firms
%
% See also PORTFOLIOS for portfolio grouping as in fama and french (1992); see also FAMAMACBETH 
%   for the cross-sectional regressions using true conditional beta; see also FF93 and VALPREM 
%   for time-series analysis of fama and french (1993)
% 
% © Lu Zhang, Inc. 2000

% initialization
lag    = 60;
[N, T] = size(R_f);

% SIMPANEL simulates 3500+ firms as in ff [93] but ff [92] only use N_f firms
% randomly draw N_f firms as a sample from the panel with 3500+ firms
% vecsub = unidrnd(N, N_f, 1);
vecsub = 1 : 1 : N_f;
R_f    = R_f(vecsub, :);               
v_f    = log(V_f(vecsub, :));            clear V_f
btm    = log(B_f(vecsub, :)) - v_f;      clear B_f vecsub

% shift the timing of size and true beta to match the convention used in simulation
btm  = [zeros(N_f, 6)    btm(:, 1 : T - 6)];

FirmIndex = zeros(N_f, T - lag - 1);        % index of portfolio sorting
Portmtrx  = zeros(100, T - lag - 1);		% matrix of portfolio returns

% some preliminary work
sizecutoff   = [772+189 236 170 144 140 128 125 119 114 124]; 
% sizecutoff = 200*ones(1, 10);

numberfirms = sum(sizecutoff);
for k = 2 : 10
	sizecutoff(k) = sizecutoff(k - 1) + sizecutoff(k);
end
sizecutoff = [1 1 + floor(sizecutoff ./numberfirms *N_f)];
sizeindex  = zeros(N_f, 1);
for k = 1 : 10
	 sizeindex(sizecutoff(k) : sizecutoff(k + 1) - 1) = 10*(k - 1);  
end

% construct a 0-1 matrix to map stocks into portfolios
M     = zeros(100, N_f);
M2    = zeros(100, N_f);
sizes = diff(sizecutoff);
for j = 0 : 99
   L        = sizes(1 + floor(j/10));		% number of stocks in a size portfolio
   stepbeta = L/10 + 1e-10;
   Imin     = (sizecutoff(1 + floor(j/10)) + mod(j,10) *stepbeta);
   Imax     = floor(Imin + stepbeta - 1);
   Imin     = floor(Imin);
   M2(j + 1, Imin:Imax) = ones(1, Imax-Imin+1);
   M (j + 1, Imin:Imax) = ones(1, Imax-Imin+1) ./(Imax-Imin+1); 
end

% sorting and exact regressions
regsize     = zeros(1, T - lag - 1);
regbtm      = zeros(1, T - lag - 1);
regbeta     = zeros(1, T - lag - 1);
regsizebtm  = zeros(2, T - lag - 1);
regsizebeta = zeros(2, T - lag - 1);
regbtmbeta  = zeros(2, T - lag - 1);

for t = lag + 2 : 1 : T
   % rebalance once every year
   if (mod(t - 1, 12) == 1) 
      % estimate prebeta first
      Y       = R_f(:, t - lag : t - 1);
	  X       = [ones(1, lag); R(t - lag : t - 1); R(t - lag - 1 : t - 2)]';   
	  regcoef = inv(X'*X)*X'*Y';
	  prebeta = (regcoef(2, :) + regcoef(3, :))';
      clear regcoef Y X;
      
      % recall size and btmvec is recorded at the beginning of a month
      sizevec = v_f(:, t);
      btmvec  = btm(:, t);
   
      % first sort by size
      Data     = [sizevec prebeta [1 : 1 : N_f]'];
      datasize = Data(:, 1);						
      [y, i]   = sort(datasize);	              clear y;		   
      Data2    = [Data(i, :) sizeindex];          clear Data;

      % Second sort by beta
      for k = 0 : 10 : 90
         J           = (Data2(:, 4) == k);
         Data3       = Data2(J, :);			% cut out a size portfolio number k
         databeta    = Data3(:, 2);	
         [y, i]      = sort(databeta);  clear y;
         Data2(J, :) = Data3(i, :);	
      end
      Data2(:, 4) = M' * [0 : 1 : 99]';
   end

   % compute portfolio returns
   FirmIndex(:, t - lag - 1) = Data2(:, 3);
   Portmtrx (:, t - lag - 1) = M * R_f(Data2(:, 3), t);

   % introduce fama-french convention: rebalancing once 12 months
   v_f(:, t) = sizevec;
   btm(:, t) = btmvec;

end

% estimate postbetas
X        = [ones(1, T - lag - 1); R(lag + 1 : T - 1); R(lag + 2 : T)]';   
regcoef  = inv(X'*X)*X'*Portmtrx';
postbeta = (regcoef(2, :) + regcoef(3, :))';
postbeta = M2'*postbeta;  clear regcoef X;

% ff regressions and cross-correlations
% initialize regression coeffients and t-stats
regffbeta     = zeros(1, T - lag - 1);	
regffsize     = zeros(1, T - lag - 1);
regffbtm      = zeros(1, T - lag - 1);
regffsizebtm  = zeros(2, T - lag - 1);
regsizeffb    = zeros(2, T - lag - 1);
regbtmffb     = zeros(2, T - lag - 1);
% initialize cross correlations
corsizeffb  = zeros(1, T - lag - 1); % size and ffbeta	
corsizebeta = zeros(1, T - lag - 1); % size and exact beta 
corbetaffb  = zeros(1, T - lag - 1); % ffbeta and exact beta
corbtmffb   = zeros(1, T - lag - 1); % ffbeta and btm
corbtmbeta  = zeros(1, T - lag - 1); % exact beta and btm
corsizebtm  = zeros(1, T - lag - 1); % size and btm

for t = 1 : T - lag - 1
   % ff regressions
   Y = R_f(:, 1 + lag + t);
   % 
   % on size alone   
   X            = [ones(N_f, 1) v_f(:, 1 + lag + t)];
   tmp          = inv(X'*X)*X'*Y;
   regffsize(t) = tmp(2);
   % on btm alone   
   X            = [ones(N_f, 1) btm(:, 1 + lag + t)];
   tmp          = inv(X'*X)*X'*Y;
   regffbtm(t)  = tmp(2);   
   % on size and btm   
   X                  = [ones(N_f, 1) v_f(: ,1 + lag + t) btm(:, 1 + lag + t)];
   tmp                = inv(X'*X)*X'*Y;
   regffsizebtm(:, t) = tmp(2 : 3);  
   % on ff-beta alone
   Y            = R_f(FirmIndex(:, t), 1 + lag + t);
   X            = [ones(N_f, 1) postbeta];
   tmp          = inv(X'*X)*X'*Y;
   regffbeta(t) = tmp(2);  
   % on size and ff-beta
   X                = [ones(N_f, 1) v_f(FirmIndex(:, t), 1 + lag + t) postbeta];
   tmp              = inv(X'*X)*X'*Y;
   regsizeffb(:, t) = tmp(2 : 3);    
   % on btm and ff-beta   
   X               = [ones(N_f, 1) btm(FirmIndex(:,t), 1 + lag + t) postbeta];
   tmp             = inv(X'*X)*X'*Y;
   regbtmffb(:, t) = tmp(2 : 3);  
   
   % calculate cross-sectional correlations
   % 
   tmp            = corrcoef(v_f(FirmIndex(:,t), 1 + lag + t), postbeta);
   corsizeffb(t)  = tmp(2);
   tmp            = corrcoef(btm(FirmIndex(:, t), 1 + lag + t), postbeta);
   corbtmffb(t)   = tmp(2);
   tmp            = corrcoef(v_f(:, 1 + lag + t), btm(:, 1 + lag + t));
   corsizebtm(t)  = tmp(2);
end; clear Data Data3 Data2;

% detail regression coefficients and t-statistics
bffbeta = mean(regffbeta) *100;
tffbeta = sqrt(T-1-lag)   *mean(regffbeta)/std(regffbeta); clear regffbeta;
bffsize = mean(regffsize) *100;
tffsize = sqrt(T-1-lag)   *mean(regffsize)/std(regffsize); clear regffsize;
bffbtm  = mean(regffbtm)  *100;
tffbtm  = sqrt(T-1-lag)   *mean(regffbtm)/std(regffbtm);   clear regffbtm;
bsizeffb     = mean(regsizeffb, 2) * 100;
tsizeffb(1)  = sqrt(T-1-lag).*mean(regsizeffb(1, :))./std(regsizeffb(1, :));
tsizeffb(2)  = sqrt(T-1-lag).*mean(regsizeffb(2, :))./std(regsizeffb(2, :)); 
clear regsizeffb;
bbtmffb      = mean(regbtmffb, 2) * 100;
tbtmffb(1)   = sqrt(T-1-lag).*mean(regbtmffb(1, :))./std(regbtmffb(1, :));
tbtmffb(2)   = sqrt(T-1-lag).*mean(regbtmffb(2, :))./std(regbtmffb(2, :));
clear regbtmffb;
bffsizebtm     = mean(regffsizebtm, 2) * 100;
tffsizebtm(1)  = sqrt(T-1-lag).*mean(regffsizebtm(1, :))./std(regffsizebtm(1, :));
tffsizebtm(2)  = sqrt(T-1-lag).*mean(regffsizebtm(2, :))./std(regffsizebtm(2, :));
clear regffsizebtm;
% detail cross-correlation
corr1 = mean(corsizeffb);        clear corsizeffb;    
corr2 = mean(corbtmffb);         clear corbtmffb;     
corr3 = mean(corsizebtm);        clear corsizebtm; 

% output data on screen and on text files
ff    = [ bffbeta   bffsize  bffbtm  bsizeffb'    bbtmffb'     bffsizebtm';
          tffbeta   tffsize  tffbtm  tsizeffb     tbtmffb      tffsizebtm  ];
corr  = [ corr1 corr2 corr3 ];
clear corr1 corr2 corr3 
