
function tbl93 = FF93(V_f, R_f, B_f, R, r, SMB, HML) 

% FF93:
%   perform on simulated panel Fama and French (1993) time series analysis:
% 
%                 R_i - r = a_i + b_i*(R - r) + s_i*SMB + h_i*HML + e_i
%
% USAGE:
%   V_f    : panel of firm size
%   R_f    : panel of firm return 
%   B_f    : panel of firm book value
%   R      : time series of market return
%   r      : time series of real interest rate
%   SMB    : small minus big
%   HML    : high minus low
% 
% OUTPUT:
%   tbl93  : 30 by 10 matrix containing the table reporting three-factor time series 
%            regressions of Fama and French (1996, table I)
%
% See FAMAMACBETH for cross-sectional regressions and VALPREM for construction of SMB and HML
%
% MODIFICATION LOG:
%   08/15/01 -- fix the bug of exogenous grouping and use the average 
% 
% © Lu Zhang, Inc. 2001

warning off

% Initializing
[N, T]  = size(R_f);
pointer = round([1; N*.20; N*.40; N*.60; N*.80; N]);  % pointers used in sorting

% shift the timing of size and book-to-market to match FF convention (and t+1 return timing)
btm     = B_f./V_f;         clear B_f;
btm     = [zeros(N, 6) btm(:, 1 : T - 6)];

% percent excess return
R_f = 100*(R_f - repmat(r, N, 1));

% further initializing
lag   = 60;
p25   = zeros([5 5 T - lag]);    % 25 portfolio returns formed on size and book-to-market

% construct 25 portfolio returns
for t = lag + 1 : T
    % sorting and rebalancing once every year
    if (mod(t - 1, 12) == 0)       
        % pick up size and btm recorded at the beginning of a month
        svec    = V_f(:, t);
        btmvec  = btm(:, t);
        IndStru = ff93Sorting(svec, btmvec, pointer);
        % construct value-weighted portfolio returns 
        tSpan            = t : t + 11;
        p25(:, :, tSpan - lag) = portConstr(IndStru, tSpan, V_f, R_f);
    end 
end

% construct Fama and French (1996) Table I
ms  = mean(p25, 3);         
sds = std(p25, 0, 3);
% Initializing
a   = zeros(5, 5);
ta  = zeros(5, 5);
b   = zeros(5, 5);
tb  = zeros(5, 5);
s   = zeros(5, 5);
ts  = zeros(5, 5);
h   = zeros(5, 5);
th  = zeros(5, 5);
R2  = zeros(5, 5);
se  = zeros(5, 5);
X   = [ones(T - lag, 1) 100*(R(lag + 1 : T) - r(lag + 1 : T))' 100*SMB' 100*HML'];   % regressors are in percent as well
% time series regressions
for i = 1 : 5
    for j = 1 : 5
        [a(i, j) ta(i, j) b(i, j) tb(i, j) s(i, j) ts(i, j) h(i, j) th(i, j) R2(i, j) se(i, j)] = ...
            olsFF93(squeeze(p25(i, j, :)), X);
    end
end
% Table I
tbl93 = [ms sds;
          a  ta;
          b  tb;
          s  ts;
          h  th;
         R2  se];  



%------------------------------subfunction declaration-----------------------------
% 
function IndStru = ff93Sorting(svec, btmvec, pointer)

% FF93SORTING:
%   given size and book-to-market vector find the indices (of stocks in the panel) for the
% 25 constructed portfolios
%
% OUTPUT:
%   IndStru : structural cell (matrix of matrices of difference lengths)
% 
% © Lu Zhang, 2001.

% Initializing
IndStru  = struct('idp11', [], 'idp12', [], 'idp13', [], 'idp14', [], 'idp15', [], ...
                  'idp21', [], 'idp22', [], 'idp23', [], 'idp24', [], 'idp25', [], ...
                  'idp31', [], 'idp32', [], 'idp33', [], 'idp34', [], 'idp35', [], ...
                  'idp41', [], 'idp42', [], 'idp43', [], 'idp44', [], 'idp45', [], ...
                  'idp51', [], 'idp52', [], 'idp53', [], 'idp54', [], 'idp55', []     );
% size sort
[ssize, sInd] = sort(svec);     clear ssize
ids1    = sInd(pointer(1)     : pointer(2));
ids2    = sInd(pointer(2) + 1 : pointer(3));
ids3    = sInd(pointer(3) + 1 : pointer(4));
ids4    = sInd(pointer(4) + 1 : pointer(5));
ids5    = sInd(pointer(5) + 1 : pointer(6));     clear sInd
% book-to-market sort
[sbtm, bInd] = sort(btmvec);    clear sbtm
idb1    = bInd(pointer(1)     : pointer(2));
idb2    = bInd(pointer(2) + 1 : pointer(3));
idb3    = bInd(pointer(3) + 1 : pointer(4));
idb4    = bInd(pointer(4) + 1 : pointer(5));
idb5    = bInd(pointer(5) + 1 : pointer(6));     clear bInd
% take index intersections to form 25 portfolios
for i = 1 : 5
    for j = 1 : 5
        eval(['IndStru.idp' num2str(i) num2str(j) '=intersect(ids' num2str(i) ', idb' num2str(j) ');']);
    end
end

%---------------------------------------------------------------------------
% 
function p25   = portConstr(IndStru, tSpan, V_f, R_f);

% PORTCONSTR:
%   construct 25 portfolios using structural cell index across the time range in tSpan
%
% OUTPUT:
%   p25 : 5 by 5 by length(tSpan) portfolio returns
% 
% © Lu Zhang, 2001.

p25 = zeros([5 5 length(tSpan)]);
for i = 1 : 5
    for j = 1 : 5
        % value-weighted returns
        eval(['p25(' num2str(i) ',' num2str(j) ', :) = sum(V_f(IndStru.idp' num2str(i) num2str(j) ...
            ', tSpan).*R_f(IndStru.idp' num2str(i) num2str(j) ', tSpan)) ./sum(V_f(IndStru.idp' ...
            num2str(i) num2str(j) ', tSpan));']);
        % equal-weighted returns
        % eval(['p25(' num2str(i) ',' num2str(j) ', :) = mean(R_f(IndStru.idp' num2str(i) num2str(j) ', tSpan));']);
    end
end

%---------------------------------------------------------------------------
% 
function [a, ta, b, tb, s, ts, h, th, R2, se] = olsFF93(Y, X)

% OLSFF93:
%   subfunction calculating single regression statistics

[T, K] = size(X);
coefs  = inv(X'*X)*X'*Y;
e      = Y - X*coefs;
R2     = 1 - var(e)/var(Y);
% adjust R2 and residual standard error by the number of regressors
R2     = 1 - (T - 1)/(T - K) *(1 - R2);
se     = sqrt(var(e)*(T - 1)/(T - K));
% tstats
covm   = e'*e /(T - K) * inv(X'*X);
tstats = coefs ./sqrt(diag(covm));
% final results
a      =  coefs(1); 
ta     = tstats(1);
b      =  coefs(2);   
tb     = tstats(2);
s      =  coefs(3);
ts     = tstats(3);
h      =  coefs(4);
th     = tstats(4);
