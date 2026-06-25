
function [sx, xinlnew] = CspSimu(xinl, xbar, rhox, stdx, T);  

% CSPSIMU:
%   Simulate an AR(1) process directly from the continuous state space
%
% USAGE:
%   xbar - long run mean
%   rhox - persistence
%   stdx - instantaneous standard deviation
%   T    - sample length
%
% © Lu Zhang, Inc. 2001

sx    = zeros(1, T + 1);
sx(1) = xinl;
innov = randn(size(sx));
xUp   = xbar + 3.5*stdx/sqrt(1 - rhox^2);
xDown = xbar - 3.5*stdx/sqrt(1 - rhox^2);

for t = 2 : T + 1
    sx(t) = min(xUp, max(xDown, xbar*(1 - rhox) + rhox*sx(t - 1) + stdx*innov(t)));
end

xinlnew = sx(end);
sx(end) = [];

% MODIFICATION LOG:
% 
%   07/27/01 -- sx(t) = min(xmax, max(xmin, xbar*(1 - rhox) + rhox*sx(t - 1) + stdx*innov(t)));
%         xmin and xmax are simply the boundary for discrete state space and cannot be used as reflection points
%         in simulations