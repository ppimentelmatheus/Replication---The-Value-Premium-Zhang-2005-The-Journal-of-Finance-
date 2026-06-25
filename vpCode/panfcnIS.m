
function [Pf, Bf, Rf, If, Af, Yf, Df, Rm, zd] = panfcn(kd0, zd0, optK, V0, k, sx, x, nx, z, N, Ts, alpha, aP, aN, delta, f, istar, rhoz, stdz)

% PANFCN:
%   bottle-neck simulation procedure in SIMPANCSP routine
%
% © Lu Zhang, Inc. 2001

kpd = kd0;
zpd = zd0;

% Initializing the firm distribution ('d' denotes values in sampling distribution)
If    = zeros(N, Ts);                                         % panel of investment
Yf    = zeros(N, Ts);                                         % panel of output
Af    = zeros(N, Ts);                                         % panel of adjustment cost
Df    = zeros(N, Ts);                                         % panel of dividend
Vf    = zeros(N, Ts);                                         % panel of firm value  
Bf    = zeros(N, Ts);                                         % panel of book value (capital stock)

for t = 1 : Ts  % be careful with the timing   
    % update current period capital --- optimal next period capital from previous period
    Bf(:, t)   = kpd;
    zd         = zpd;

    % update the cross-section of next period capital stock and firm value
    % 
    % step 1 -- linear interpolation (extrapolation) along x dimension
    index  = find(x >= sx(t));
    if length(index) == 0, 
        % sx(t) is larger than all points on x grid
        optKix = squeeze(optK(:, end, :) + ((sx(t) - x(end))/(x(end) - x(end - 1)))*(optK(:, end, :) - optK(:, end - 1, :)));
        V0ix   = squeeze(V0(:, end, :) + ((sx(t) - x(end))/(x(end) - x(end - 1)))*(V0(:, end, :) - V0(:, end - 1, :)));
    elseif length(index) == nx
        % sx(t) is smaller than all points on x grid 
        optKix = squeeze(optK(:, 1, :) - ((x(1) - sx(t))/(x(2) - x(1)))*(optK(:, 2, :) - optK(:, 1, :)));
        V0ix   = squeeze(V0(:, 1, :) - ((x(1) - sx(t))/(x(2) - x(1)))*(V0(:, 2, :) - V0(:, 1, :)));
    else
        % sx(t) lies in the middle of x grid    
        down   = index(1);
        up     = down - 1;
        pup    = (x(down) - sx(t))/(x(down) - x(up));
        optKix = squeeze(pup*optK(:, up, :) + (1 - pup)*optK(:, down, :));
        V0ix   = squeeze(pup*V0(:, up, :) + (1 - pup)*V0(:, down, :));
    end
    % 
    % step 2 -- point-by-point bilinear interpolation (extrapolation) along kd and zd dimension
    for j = 1 : N
        % find up and down index of Bf(j, t) on k grid
        index = find(k >= Bf(j, t));
        down  = index(1);
        if down == 1, % kd(j) == kmin
            optKixk = optKix(1, :);
            V0ixk   = V0ix(1, :);
        else
            % linear interpolation along k dimension (continuous space)
            up      = down - 1;
            pup     = (k(down) - Bf(j, t))/(k(down) - k(up));
            optKixk = pup*optKix(up, :) + (1 - pup)*optKix(down, :);
            V0ixk   = pup*V0ix(up, :) + (1 - pup)*V0ix(down, :);
        end
        % linear interpolation along z dimension (continuous space)
        kpd(j)   = interp1(z, optKixk', zd(j), 'linear', 'extrap');
        Vf(j, t) = interp1(z, V0ixk', zd(j), 'linear', 'extrap');
    end
    
    % update cross-sectional dividends
    If(:, t) = kpd - (1 - delta)*Bf(:, t);  
    tmpIf    = If(:, t)./Bf(:, t) - istar;
    Af(:, t) = ((aP*(tmpIf >= 0) + aN*(tmpIf < 0))/2).*(tmpIf.^2).*Bf(:, t);  
    Yf(:, t) = exp(sx(t) + zd).*(Bf(:, t).^alpha);
    
    % update idiosyncratic shock in continuous state space
    shock  = randn(N/2, 1);
    zpd    = max(-3.5*(stdz/sqrt(1 - rhoz^2)), min(3.5*(stdz/sqrt(1 - rhoz^2)), rhoz*zd + stdz*[shock; -shock]));
end

% dividend
Df  = Yf - If - Af - f;
% important step: redefined firm value to be ex dividend to conform to empirical studies
Pf  = Vf - Df;      
% cross-sectional stock return
Rf  = (Pf(:, 2:end) + Df(:, 2:end)) ./Pf(:, 1:end-1);     % ex dividend convention
% value-weighted market return
Rm  = sum(Pf(:, 1:end-1).*Rf) ./sum(Pf(:, 1:end-1));     
