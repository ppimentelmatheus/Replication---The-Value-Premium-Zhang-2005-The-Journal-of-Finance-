% PLOTPOLICY:
%   Plot policy function for presentations
% 
% © Lu Zhang, Inc. 2001

clear all; clc; format compact; format short

load plot_data
load Params

nx = length(x);
nz = length(z);
nh = length(h);
nk = length(k);

% reshape to facilitate 3-D plotting
optK = reshape(optK, [nk nh nx nz]);
V    = reshape(V, [nk nh nx nz]);
I    = reshape(I, [nk nh nx nz]);
div  = reshape(div, [nk nh nx nz]);

% imx  = find(x == mean(x));
% imz  = find(z == 0);
% imk  = find(k == 1);
% imh  = find(h == mean(h));

% value function plots
figure(1); 
mesh(k, h, squeeze(V(:, :, imx, imz))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('P', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_kh.eps
figure(2); 
mesh(k, x, squeeze(V(:, imh, :, imz))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('x', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_kx.eps
figure(3); 
mesh(k, z, squeeze(V(:, imh, imx, :))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_kz.eps
figure(4); 
mesh(h, x, squeeze(V(imk, :, :, imz))', 'LineWidth', 2)
xlabel('P', 'FontS', 15); ylabel('x', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_hx.eps
figure(5); 
mesh(h, z, squeeze(V(imk, :, imx, :))', 'LineWidth', 2)
xlabel('P', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_hz.eps
figure(6); 
mesh(x, z, squeeze(V(imk, imh, :, :))', 'LineWidth', 2)
xlabel('x', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('V', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\v_xz.eps

% optimal investment plots
figure(11); 
mesh(k, h, squeeze(I(:, :, imx, imz))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('P', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_kh.eps
figure(12); 
mesh(k, x, squeeze(I(:, imh, :, imz))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('x', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_kx.eps
figure(13); 
mesh(k, z, squeeze(I(:, imh, imx, :))', 'LineWidth', 2)
xlabel('K', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_kz.eps
figure(14); 
mesh(h, x, squeeze(I(imk, :, :, imz))', 'LineWidth', 2)
xlabel('P', 'FontS', 15); ylabel('x', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_hx.eps
figure(15); 
mesh(h, z, squeeze(I(imk, :, imx, :))', 'LineWidth', 2)
xlabel('P', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_hz.eps
figure(16); 
mesh(x, z, squeeze(I(imk, imh, :, :))', 'LineWidth', 2)
xlabel('x', 'FontS', 15); ylabel('z', 'FontS', 15); zlabel('I', 'FontS', 15)
set(gca, 'FontS', 15);
print -deps e:\Research\MyPapers\ValPrem\Documents\FIGURES\i_xz.eps
