
function [tbl99, pSMB, pHML, pSH, pBH, pSL, pBL] = ValPrem(V_f, R_f, B_f, R, r) 

% VALPREM: Perform on simulated panel:
%     (1) construction of SMB and HML
%     (2) summary statistics for monthly percent three-factor explanatory returns
%
% USAGE:
%   V_f    : panel of firm size
%   R_f    : panel of firm return
%   B_f    : panel of firm book value
%   R      : time series of market return
%   r      : time series of real interest rate
% 
% OUTPUT:
%   pSMB   : time series of SMB portfolio returns
%   pHML   : time series of HML portfolio returns
%   tbl99  : table reporting summary statistics of SMB and HML with format:
% 
%                       excess return    SMP    HML   %   SM    SH    BL    BM    BH
%       percent mean:
%       percent  std:
%       t-statistics:
%
% See FAMAMACBETH for cross-sectional regressions and FF93 for time series regressions
% 
% © Lu Zhang, 2001.

% Initializing
[N, T] = size(R_f);

% shift the timing of size and book-to-market to match FF convention (and t+1 return timing)
btm  = B_f./V_f;       clear B_f;
btm  = [zeros(N, 6)    btm(:, 1 : T - 6)];

% percentage of firms in size and btm portfolios -- capturing ff 93 practice
small  = round(0.50*N);
low    = round(0.30*N);
lowmed = round(0.70*N);

% further initializing
lag = 60;
pSL = zeros(1, T - lag);     % portfolio returns
pSM = zeros(1, T - lag);
pSH = zeros(1, T - lag);
pBL = zeros(1, T - lag);
pBM = zeros(1, T - lag);
pBH = zeros(1, T - lag);

% construct portfolio returns
for t = lag + 1 : T
    % sorting and rebalancing once every year
    if (mod(t - 1, 12) == 0)       
        % pick up size and btm recorded at the beginning of a month
        svec    = V_f(:, t);
        btmvec  = btm(:, t);

        if 1    % ff 93 cut
           % size split
           [ssize, index] = sort(svec);      clear ssize
           idS     = index(1 : small);
           idB     = index(small + 1 : end);
           % book-to-market sort
           [sbtm, index] = sort(btmvec);     clear sbtm
           idL     = index(1 : low);
           idM     = index(low + 1 : lowmed);
           idH     = index(lowmed + 1 : end);
        else    % endogenous cut
           % size split
           bpSize  = median(svec);          
           idS     = find(svec <= bpSize);           % index of small stocks
           idB     = find(svec >  bpSize);           % index of large stocks
           % book-to-market sort
           bpLow   = min(btmvec) + (max(btmvec) - min(btmvec))*0.30;
           bpMed   = min(btmvec) + (max(btmvec) - min(btmvec))*0.70;
           idL     = find(btmvec <= bpLow);
           idM     = find(btmvec >  bpLow & btmvec <= bpMed);
           idH     = find(btmvec >  bpMed);
        end
        
        % take index intersections to form six portfolios
        idSL    = intersect(idS, idL);
        idSM    = intersect(idS, idM);
        idSH    = intersect(idS, idH);
        idBL    = intersect(idB, idL);
        idBM    = intersect(idB, idM);
        idBH    = intersect(idB, idH);
        % construct value-weighted portfolio returns
        tSpan   = t - lag : t - lag + 11;
        pSL(tSpan) = sum(V_f(idSL, tSpan).*R_f(idSL, tSpan)) ./sum(V_f(idSL, tSpan));
        pSM(tSpan) = sum(V_f(idSM, tSpan).*R_f(idSM, tSpan)) ./sum(V_f(idSM, tSpan));
        pSH(tSpan) = sum(V_f(idSH, tSpan).*R_f(idSH, tSpan)) ./sum(V_f(idSH, tSpan));
        pBL(tSpan) = sum(V_f(idBL, tSpan).*R_f(idBL, tSpan)) ./sum(V_f(idBL, tSpan));
        pBM(tSpan) = sum(V_f(idBM, tSpan).*R_f(idBM, tSpan)) ./sum(V_f(idBM, tSpan));
        pBH(tSpan) = sum(V_f(idBH, tSpan).*R_f(idBH, tSpan)) ./sum(V_f(idBH, tSpan));
    end
end


% SMB and HML
pSMB        = (pSL + pSM + pSH)/3 - (pBL + pBM + pBH)/3;
pHML        = (pSH + pBH)/2 - (pSL + pBL)/2;

% some visual
if 0
meanS       = [mean(pSL)^12   mean(pSM)^12   mean(pSH)^12]
meanB       = [mean(pBL)^12   mean(pBM)^12   mean(pBH)^12]
meanH       = [mean(pSH)^12   mean(pBH)^12]
meanL       = [mean(pSL)^12   mean(pBL)^12]
end

% construct Fama and French (1999) Table I
tbl99       = zeros(3, 9);
% tbl99     = zeros(3, 3);
tbl99(1, :) = 100*[ mean(R - r)   mean(pSMB)   mean(pHML) ...
                       mean(pSL-1)   mean(pSM-1)  mean(pSH-1) mean(pBL-1) mean(pBM-1) mean(pBH-1) ];
tbl99(2, :) = 100*[ std(R - r)  std(pSMB)  std(pHML) ...
                             std(pSL)    std(pSM)   std(pSH)  std(pBL)  std(pBM)  std(pBH) ];
tbl99(3, :) = sqrt(T) * tbl99(1, :)./ tbl99(2, :);    % t-statistics            


