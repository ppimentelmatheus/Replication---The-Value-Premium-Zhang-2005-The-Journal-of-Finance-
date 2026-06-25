
function btmPanel = ponShall(Rm, btm)

% PONSHALL: 
%   perform time series regression of market return on aggregate book-to-market
% as in Pontiff and Schall (1999)
% 
% USAGE:
%   Rm  -- 1 by T (multiply of 840) simulated time series of market return 
%   btm -- 1 by T (multiply of 840) simulated time series of book-to-market
% 
% NOTES: Rm(t) and btm(t) are already lined up in that btm(t) is measured at the
%   beginning of period t and Rm(t) is from t to t+1. Although I always defined
%   Rm(t+1) to be the return from t to t+1 in theory but in practice the first 
%   column of P gets lost so the first Rm(1) is computed using (P(2) + D(2))/P(1)!
% 
% © Lu Zhang, Inc. 2001.


T    = length(Rm);
plen = 840;
nsim = T/plen;
if mod(T, plen) ~= 0
    disp('The length of Rm and btm is not right!')
    btmPanMean = [];
    btmPanStd  = [];
    return
end
Rm  = reshape(Rm, plen, nsim)';
btm = reshape(btm, plen, nsim)';

% perform Psim ontiff and Schall [1998] time series analysis
btmPanel = zeros(2, 5);                        % result matrix
col      = 1 : 12 : plen;                      % from monthly btm to annual btm
for is = 1 : nsim
    % monthly frequency
    btmMonth   = log(btm(is, :))';
    retMonth   = Rm(is, :)';
    % [bm, tm, Rsqm] = linregAdjust(retMonth, [ones(plen, 1) btmMonth]);  
    [bm, tm, Rsqm, chi2, p] = MultiRegressNw(retMonth, [ones(plen, 1) btmMonth], [], [], 0);
    
    % annual frequency
    btmAnnual = btmMonth(col);
    retAnnual = reshape(retMonth, 12, plen/12); 
    retAnnual = prod(retAnnual, 1)';
    % [ba, ta, Rsqa] = linregAdjust(retAnnual, [ones(plen/12, 1) btmAnnual]);  
    [ba, ta, Rsqa, chi2, p] = MultiRegressNw(retAnnual, [ones(plen/12, 1) btmAnnual], [], [], 0);
   
    % store simulated statistics on dividend yield and book-to-market
    btmPanel = [ bm' tm' Rsqm;
                 ba' ta' Rsqa ];
end

% report results
% btmPanMean = mean(btmPanel, 3);
% btmPanStd  = std(btmPanel, 0, 3);