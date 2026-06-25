
function [Junk] = plotFF95(roevec, btmvec, ROElow, ROEhigh)

% PLOTFF95
%   plot figures in FF 95 using the results computed in FF95 routine
%
% © Lu Zhang, Inc. 2001

Junk = [];
xseries = -5 : 1: 5;

% plot profitability: Figure 1 in FF 95
figure(951)
plot(xseries, roevec(:, 1), '-', 'LineWidth', 2); hold on; grid on;
plot(xseries, roevec(:, 2), ':', 'LineWidth', 2); hold off;
xlabel('Formation Year', 'FontS', 15); ylabel('Profitability', 'FontS', 15);
set(gca, 'FontS', 15);
legend('Growth', 'Value'); axis([-6 6 -0.05 0.25])
% print -deps e:\Research\ValPrem\Documents\FIGURES\roe11.eps

% plot book-to-market: Figure 2 in FF 95
figure(952)
plot(xseries, btmvec(:, 1), '-', 'LineWidth', 2); hold on; grid on;
plot(xseries, btmvec(:, 2), ':', 'LineWidth', 2); hold off;
xlabel('Formation Year', 'FontS', 15); ylabel('BE/ME', 'FontS', 15);
set(gca, 'FontS', 15);
legend('Growth', 'Value', 0); axis([-6 6 0.50 1.10])
% print -deps e:\Research\ValPrem\Documents\FIGURES\btm11.eps

return

% plot profitability: Figure 3 in FF 95
figure(953);
plot(ROElow(2:end), '-', 'LineWidth', 2); hold on; grid on
plot(ROEhigh(2:end), ':', 'LineWidth', 2); hold off; 
xlabel('Time Series', 'FontS', 15); ylabel('Profitability', 'FontS', 15);
set(gca, 'FontS', 15);
legend('Low B/M', 'High B/M');
print -deps e:\Research\ValPrem\Documents\FIGURES\profitability.eps

