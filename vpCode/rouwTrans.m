
function [transP, dscSp] = rouwTrans(rho, ave, dev, num)

% ROUWTRANS:
%    Construct transition probability matrix for discretizing an AR(1) process
% This procedure is from Rouwenhorst (1995) which can calibrate an extremely
% persistent AR(1) process that Tauchen and Hussey (1991) procedure cannot handle
% 
% USAGE:
%   rho - persistent level of the process (close to one)
%   ave - mean serving as the middle point of the discrete state space
%   dev - step size of the even-spaced grid 
%   num - num of grid point on the discretized process
% 
% OUTPUT:
%   dscSp  - discrete state space (num by 1 vector)
%   transP - transition probability matrix over the grid with sum(transP) == 1
%        
% NOTES: 
%   (1) How to pick rho and dev to match an AR(1) with persistence rho and conditional
%       volatility sigma? First, rho is itself. That's easy. Now consider the choice of dev.
%       Rouwenhorst (1995) shows that the unconditional volatility of discrete process is: 
%         
%                                epsilon/sqrt(num - 1) 
% 
%       where epsilon is from [ave - epsilon ave + epsilon] the grid
%       interval (now you see why num has to be an odd number so as to keep the interval symmetric 
%       around ave!). Next, using the following two facts: (i) epsilon = (num - 1)/2 *dev and (ii)
%       conditional volatility = unconditional volatility * sqrt(1 - rho^2), I have the cook recipe:
% 
%               dev = 2*sigma^{cond} /sqrt((1 - rho^2)*(num - 1))
%       
%   (2) To go from quarterly rho_q and sigma_q (conditional volatility) to monthly counterparts,
%       simply do:      rho_m = rho_q^(1/4)    and sigma_m = sigma_q/4
%       where the latter equality follows since conditional shocks are IID by definition
% 
% © Lu Zhang, Inc. 2001

% Quality Control
if mod(num, 2) == 0,
    fprintf('ROUWTRANS only support odd number of grid points: choose num again!');
    return
end

% Discrete state space
dscSp = linspace(ave-(num-1)/2*dev, ave+(num-1)/2*dev, num)';

% Transition probability matrix
p = (rho + 1)/2;
q = p;
transP = [ p^2       p*(1-q)         (1-q)^2  ;
           2*p*(1-p) p*q+(1-p)*(1-q) 2*q*(1-q);
           (1-p)^2   (1-p)*q         q^2        ]';
       
while size(transP, 1) <= num - 1
    len    = size(transP, 1); 
    % Rule 1: See Rouwenhorst (1995)
    transP = p*[ transP zeros(len, 1); zeros(1, len) 0 ] + ...
                 (1 - p)*[zeros(len, 1) transP; 0 zeros(1, len)] + ...
                 (1 - q)*[zeros(1, len) 0; transP zeros(len, 1)] + ...
                 q*[0 zeros(1, len); zeros(len, 1) transP];
    % Rule 2: See Rouwenhorst (1995)
    transP(2:end-1, :) = transP(2:end-1, :)/2;
end

% Take transpose to convert transP to my convention
transP = transP';

% Some Quality Control
if max(abs(sum(transP) - ones(1, num))) >= 1e-8,
    max(abs(sum(transP) - ones(1, num)))
    warning('Problem in ROUWTRANS routine!');
end
