function robustlinreg_demo
%ROBUSTLINREG_DEMO Demo of robust linear regression
%
%   robustlinreg_demo;
%
%       This is a small program that demonstrates the use of robust
%       linear regression.
%

% Created by Dahua Lin, on Jan 15, 2012
%

%% prepare data

n = 500;
noise = 0.05;
outlier_fraction = 0.3;
outliers = randpick(n, round(n * outlier_fraction));

a0 = rand(2, 1);

x = randn(1, n);
Xa = [x; ones(1, n)];
y = a0' * Xa + randn(1, n) * noise;

out_dev = 5;
y(outliers) = y(outliers) + rand(1, numel(outliers)) * out_dev; 


%% regression
 
% generalized linear regression using quadratic loss
f1 = genlinregf(x, y, [], @quadloss, 1e-4);
a1 = bfgsfmin(f1, zeros(2,1));

% generalized linear regression using huber loss
f2 = genlinregf(x, y, [], @(e) huberloss(e, 1e-2), 1e-4);
a2 = bfgsfmin(f2, a1);

% generalized linear regression using bisquare loss
f3 = genlinregf(x, y, [], @(e) bisquareloss(e, 1), 1e-4);
a3 = bfgsfmin(f3, a2);

write_line('true', a0);
write_line('quad loss', a1); 
write_line('huber loss', a2);
write_line('bisquare loss', a3);


%% visualization

xmin = min(x);
xmax = max(x);

figure;
plot(x, y, '.');
hold on; plot_line(a1, xmin, xmax, 'r-');
hold on; plot_line(a2, xmin, xmax, 'm-');
hold on; plot_line(a3, xmin, xmax, 'g-');


%% auxiliary functions

function plot_line(a, xmin, xmax, symb)

x0 = xmin;
x1 = xmax;
y0 = a(1) * x0 + a(2);
y1 = a(1) * x1 + a(2);

plot([x0 x1], [y0 y1], symb, 'LineWidth', 1.5);

function write_line(name, a)

fprintf('line (%s): y = %.4f x + %.4f\n', name, a(1), a(2));





