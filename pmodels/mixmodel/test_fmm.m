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

function test_fmm_on_gaussgm(K, n)

% model config

d = 2;
gm = gaussgm(d);
gm.cov_type = 'full';
gm.tie_cov = false;

mu = 5 * randn(d, K);
C0 = zeros(d, d, K);
for k = 1 : K   
    R = orth(randn(d,d));
    cc = R * diag([1 0.3]) * R';
    C0(:,:,k) = (cc + cc') * 0.5;    
end

% data

X = cell(1, K);
for k = 1 : K
    cg = gaussd.from_mp(mu(:,k), gsymat(C0(:,:,k)));
    X{k} = cg.sample(n);
end
X = [X{:}];

% estimate

FM = fmm_em(gm, X, K);

% visualize

figure;
plot(X(1,:), X(2,:), 'b+', 'MarkerSize', 5);

models = FM.params;
hold on;
plot_ellipse(models, 1, 'r-', 'LineWidth', 2);
plot_ellipse(models, 3, 'r-', 'LineWidth', 1);

axis equal;



function test_fmm_on_gaussgm_gp(K, n)

% model config

d = 2;
prior = gaussd.from_mp('s', zeros(d, 1), 25, 'ip');
gm = gaussgm_gp(prior, 1);

% data

thetas = gm.prior.sample(K);
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

