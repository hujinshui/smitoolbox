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

a0 = rand(2, 1);, 

x = rand(n, 1);
X = [x, ones(n, 1)];
y = X * a0 + randn(n, 1) * noise;

sout = randpick(n, round(n * outr));
y(sout) = y(sout) + rand(numel(sout), 1) * 5; 


figure;
plot(x, y, '.');
hold on; plot_line(a0, 'r-');

a = llsq(X, y);
hold on; plot_line(a, 'g--');


disp('Solve using fminunc:');
f = rbstlr(X, y, [], 'bisquare', [], 0.01);
options = optimset('LargeScale', 'on', 'GradObj', 'on', ...
    'TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'iter');
ar1 = fminunc(f, zeros(2,1), options);
hold on; plot_line(ar1, 'm-');


function plot_line(a, symb)

plot([0 1], [a(2), a(1)+a(2)], symb, 'LineWidth', 1.5);