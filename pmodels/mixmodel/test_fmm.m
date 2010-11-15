function test_fmm(model_type, K, n)
% A simple script to demo how to use finite mixture model estimation
%
%   test_fmm(model_type, K, n);
%
%   Here, model_type can be either of the following values:
%   - 'gaussgm'
%   - 'gaussgm_gp'
%
%   K:  the number of component models
%   n:  the number of sample from each component
%

%% main

switch model_type
    case 'gaussgm'
        test_fmm_on_gaussgm(K, n);
    case 'gaussgm_gp'
        test_fmm_on_gaussgm_gp(K, n);
    otherwise
        error('test_fmm:invalidarg', 'The model_type is invalid.');
end


%% core functions

function test_fmm_on_gaussgm_gp(K, n)

% model config

d = 2;
prior = gaussd.from_mp(0, udmat(d, 25), 'ip');
noise_J = udmat(d, 1);
gm = gaussgm_gp(prior, [], noise_J);

% data

thetas = gm.pri_sample(K);
X = cell(1, K);
for k = 1 : K
    X{k} = gm.sample(thetas(:,k), n);
end
X = [X{:}];

% estimate

FM = fmm_em(gm, X, K);

% visualize

figure;
plot(X(1,:), X(2,:), 'b+', 'MarkerSize', 5);

models = gm.param_models(FM.params, 'mp');
hold on;
plot_ellipse(models, 1, 'r-', 'LineWidth', 2);
plot_ellipse(models, 3, 'r-', 'LineWidth', 1);

axis equal;
