function gaussgm_demo()
% A script to demo the inference over Gaussian generative model
%
%   gaussgm_demo;
%

% Created by Dahua Lin, on Nov 14, 2010
%


%% model configuration

prior = gaussd.from_mp('f', zeros(2, 1), eye(2) * 1e4, 'ip');

sigma = 0.1;
gm = gaussgm_gp(prior, sigma);

%% generate data

n = 60;
theta = randn(2, 1);
x = bsxfun(@plus, theta, gm.gnoise.sample(n));

%% do inference

theta_map = gm.estimate_map(x);
sp = gm.pos_sample(x, 1, 50);


%% visualize

figure;

plot(x(1,:), x(2,:), 'b.', 'MarkerSize', 12);

hold on;
plot(theta_map(1), theta_map(2), 'ro', 'MarkerSize', 20, 'LineWidth', 2);

hold on;
plot(sp(1,:), sp(2,:), 'm+', 'MarkerSize', 10);

axis equal;
