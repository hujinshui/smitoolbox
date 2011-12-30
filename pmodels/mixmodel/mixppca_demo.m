function mixppca_demo(n, K)
% A program to demonstrate the use of mixture of PPCA
%
%   mixppca_demo(n, K);
%
%       Here, n is the number of sample points, K is the number of 
%       mixture components to be used
%
%       By default n = 50000, and K = 8;
%

% Created by Dahua Lin, on Nov 6, 2011
%

%% get input

if nargin < 1
    n = 50000;
end

if nargin < 2
    K = 8;
end


%% generate data

sigma0 = 0.1;
r0 = 2;
rgn = [-3, 3, -3, 3];

theta = rand(1, n) * (2 * pi);
r = r0 + randn(1, n) * sigma0;

X = [r .* cos(theta); r .* sin(theta)];

%% Do estimation

gm = ppca_gm(2, 1);
c0 = 1;
state = fmm_std('em', gm, [], c0);
Z0 = ceil(theta / (2 * pi) * K);
state = state.initialize_by_group(X, [], K, Z0);

opts = varinfer_options([], ...
    'maxiters', 50, ...
    'display', 'eval', ...
    'tol', 1e-6);

R = varinfer_drive(state, opts);
models = R.sol.params;

%% visualize

% models

figure;
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 3);

axis(rgn);
axis equal;

for k = 1 : K
    gk = gaussd('m', models(k));
    hold on;
    gaussd_ellipse(gk, 2, 500, 'r', 'LineWidth', 2);
end

