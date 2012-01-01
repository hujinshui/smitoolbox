function ica_demo(method)
% A simple program to demonstrate how ICA is used to decouple noisy signals
%
%   ica_demo(method);
%
%   You can choose the method of running ICA: 'defl' or 'symm'.
%

% Created by Dahua Lin, on Dec 31, 2011


%% prepare signals

t = (0 : 0.02 : 2 * pi).';
d = length(t);
n = 500;

X1 = bsxfun(@plus, sin(t), randn(d, n) * 1e-1);
X2 = bsxfun(@plus, cos(t), randn(d, n) * 1e-1);

X = [X1 X2];
sX = sqrt(mean(sum(X.^2, 1)));

%% run Fast ICA

W = ica_fast(X, 2, [], 'method', method, 'verbose', true);

%% visualize

figure;
plot(t, X, 'g.', 'MarkerSize', 3);
hold on;
plot(t, W * sX, 'r-', 'LineWidth', 2);


