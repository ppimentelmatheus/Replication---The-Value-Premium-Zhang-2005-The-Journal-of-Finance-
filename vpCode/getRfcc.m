
function rf = getRfcc(x, Qx, beta, eta_1, eta_2)

% GETRFCC:
%   Obtain Closed-form Real Interest Rate Using Discrete State Space 
% Transition Matrix and Parametric Pricing Kernel Process:
% 
%                         log M_{t,t+1} = log(beta) + eta (x_t - x_{t+1})
%   and
%                                  eta = eta_1 + eta_2 (x_t - xbar)
% 
%   where eta_1 > 0 and eta_2 < 0 implying risk-aversion rises in recession and shrinks in expansions
%
% © Lu Zhang, 2001.

xbar  = mean(x);
nx    = length(x);
% Construct EM
EM    = zeros(1, nx); 
for ix = 1 : nx
    EM(ix) = beta*exp((eta_1 - eta_2*xbar)*x(ix) + eta_2*x(ix)^2)*...
        (exp(-(eta_1 - eta_2*xbar + eta_2*x(ix))*x')*Qx(:, ix));
end
% Real Interest Rate
rf    = 1./EM;