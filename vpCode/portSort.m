
function [p25] = portSort(Vf, Bf, Rf, r)

% PORTSORT: Construct Fama and French 25 portfolios, 10 B/M portfolios, and 10 Size portfolios
% 
% INPUT:
%   Vf: N by T matrix of firm capitalization
%   Bf: N by T matrix of book value of firms
%   Rf: N by T matrix of individual firm return
% 
% OUTPUT:
%   p25: 
% 
% © Lu Zhang, 2001.

% warning off

% Initializing
[N, T]  = size(Rf);
pointer = round([1; N*.20; N*.40; N*.60; N*.80; N]);  % pointers used in sorting

% shift the timing of size and book-to-market to match FF convention (and t+1 return timing)
btm     = Bf./Vf;         clear Bf;
btm     = [zeros(N, 6) btm(:, 1 : T - 6)];

% percent excess return of stocks
Rf = 100*(Rf - repmat(r, N, 1));

% further initializing
lag   = 0;
p25   = zeros([5 5 T - lag]);    % 25 portfolio returns formed on size and book-to-market

% construct 25 portfolio returns
for t = lag + 1 : T
    % sorting and rebalancing once every year
    if (mod(t - 1, 12) == 0)       
        % pick up size and btm recorded at the beginning of a month
        svec    = Vf(:, t);
        btmvec  = btm(:, t);
        IndStru = ff93Sorting(svec, btmvec, pointer);
        % construct value-weighted ff25 portfolio returns 
        tSpan            = t : t + 11;
        p25(:, :, tSpan - lag) = p25Constr(IndStru, tSpan, Vf, Rf);
    end 
end


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
function p25   = p25Constr(IndStru, tSpan, V_f, R_f);

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
