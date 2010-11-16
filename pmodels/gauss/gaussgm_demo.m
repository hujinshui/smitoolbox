function gaussgm_demo()
% A script to demo the inference over Gaussian generative model
%
%   gaussgm_demo;
%

% Created by Dahua Lin, on Nov 14, 2010
%


%% model configuration

A = randn(2, 2);
b = randn(2, 1);

prior = gaussd.from_mp(0, gsymat(eye(2) * 1e+6), 'ip');

sigma = udmat(2, 0.1);
noise = gaussd.from_mp(0, sigma);

gm = gaussgm_gp(prior, A, inv(sigma));

%% generate data

n = 100;
theta = randn(2, 1);
x = bsxfun(@plus, A * theta + b, noise.sample(n));

%% do inference

y = bsxfun(@minus, x, b);

pos = gm.get_posterior(y).to_gaussd('mp');
emap = gm.estimate_map(y);

sp = gm.pos_sample(y, 1, [], 50);


%% visualize

figure;
% plot(x(1,:), x(2,:), 'b+', 'MarkerSize', 5);

hold on;
plot(sp(1,:), sp(2,:), 'm+', 'MarkerSize', 5);

hold on;
plot(emap(1,:), emap(2,:), 'b+', 'MarkerSize', 10);

hold on;
plot(theta(1), theta(2), 'r.', 'MarkerSize', 20);

plot_ellipse(pos, 1, 'm-');
plot_ellipse(pos, 3, 'm-');

axis equal;
