function rbstlr_demo()
% A simple script demonstrating the use of rbstlr
%
%   rbstlr_demo;
%

% Created by Dahua Lin, on Jan 6, 2011
%

n = 500;
noise = 0.05;
outr = 0.2;

a0 = rand(2, 1);

x = rand(n, 1);
X = [x, ones(n, 1)];
y = X * a0 + randn(n, 1) * noise;

sout = randpick(n, round(n * outr));
y(sout) = y(sout) + rand(numel(sout), 1) * 5; 


% figure;
plot(x, y, '.');
hold on; plot_line(a0, 'r-');

a = llsq(X, y);
hold on; plot_line(a, 'g--');

disp('Solve using Newton:');
ar1 = rbstlr(X, y, 'bisquare', [], ...
    'Method', 'newton', 'Display', 'iter', 'L2R', 0.01);
disp(' ');

disp('Solve using BFGS:');
ar2 = rbstlr(X, y, 'bisquare', [], ...
    'Method', 'bfgs', 'Display', 'iter', 'L2R', 0.01);
disp(' ');

disp('Solve using IRLS:');
ar3 = rbstlr(X, y, 'bisquare', [], ...
    'Method', 'irls', 'Display', 'iter', 'L2R', 0.01);
disp(' ');


hold on; plot_line(ar1, 'm-');

fprintf('ar (Newton) = %.4f x + %.4f\n', ar1(1), ar1(2));
fprintf('ar (BFGS)   = %.4f x + %.4f, diff = %.3g\n', ...
    ar2(1), ar2(2), Linfdiff(ar1, ar2));
fprintf('ar (IRLS)   = %.4f x + %.4f, diff = %.3g\n', ...
    ar3(1), ar3(2), Linfdiff(ar1, ar3));
disp(' ');


function plot_line(a, symb)

plot([0 1], [a(2), a(1)+a(2)], symb, 'LineWidth', 1.5);