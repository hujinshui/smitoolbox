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

a0 = rand();
b0 = rand();

x = randn(1, n);
y = (a0 * x + b0) + randn(1, n) * noise;

out_dev = 5;
y(outliers) = y(outliers) + rand(1, numel(outliers)) * out_dev; 


%% regression
 
% generalized linear regression using quadratic loss
[a1, b1] = glinreg(x, y, [], @quadloss, 1e-4);

% generalized linear regression using huber loss
[a2, b2] = glinreg(x, y, [], @(e) huberloss(e, 1e-2), 1e-4);

% generalized linear regression using bisquare loss
[a3, b3] = glinreg(x, y, [], @(e) bisquareloss(e, 1), 1e-4);

write_line('true', a0, b0);
write_line('quad loss', a1, b1); 
write_line('huber loss', a2, b2);
write_line('bisquare loss', a3, b3);


%% visualization

xmin = min(x);
xmax = max(x);

figure;
plot(x, y, '.');
hold on; plot_line(a1, b1, xmin, xmax, 'r-');
hold on; plot_line(a2, b2, xmin, xmax, 'm-');
hold on; plot_line(a3, b3, xmin, xmax, 'g-');


%% auxiliary functions

function plot_line(a, b, xmin, xmax, symb)

x0 = xmin;
x1 = xmax;
y0 = a * x0 + b;
y1 = a * x1 + b;

plot([x0 x1], [y0 y1], symb, 'LineWidth', 1.5);

function write_line(name, a, b)

fprintf('line (%s): y = %.4f x + %.4f\n', name, a, b);





