function R = glmm_demo(method, K)
% A script to demo the inference over mixture of Gaussian linear models
%
%   R = glmm_demo(method, K);
%   
%       Input arguments:
%       - method:       either of the following strings
%                       - 'em': expectation-maximization
%                       - 'hard-em': expectation-maximization (hard)
%                       - 'gibbs': Gibbs sampling
%
%       - K:            the number of mixture components
%

%   History
%   -------
%       - Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~ischar(method)
    error('glmm_demo:invalidarg', 'The method should be a string.');
end
method = lower(method);

if ~(isnumeric(K) && isscalar(K) && K >= 2)
    error('glmm_demo:invalidarg', ...
        'K should be a positive integer with K >= 2.');
end

%% prepare data

n = 1000;   % # samples / class

pri_sigma = 10 * sqrt(K);
g0 = gaussd('m', 0, pdmat('s', 2, pri_sigma));

Cx = pdmat('s', 2, 1);
U0 = gaussd_sample(g0, K);

Xs = cell(1, K);
for k = 1 : K
    Xs{k} = gaussd_sample(gaussd('m', U0(:,k), Cx), n);    
end
X = [Xs{:}];

%% fit models

Jx = pdmat_inv(Cx);
glm = gauss_lingen(Jx);
gpri = gausspri(g0);

pri_count = 1;
state = fmm_std(method, glm, gpri, pri_count);
L0 = randi(K, 1, size(X,2));
state = state.initialize_by_group(X, [], K, L0);

switch method
    case {'em', 'hard-em'}
        opts = varinfer_options([], ...
            'maxiters', 200, 'tol', 1e-6, 'display', 'eval');
        R = varinfer_drive(state, opts);
        
    otherwise
        error('The method %s is not supported', method);
end

%% Visualize

U = R.sol.params;

switch method
    case 'em'
        [~, Zm] = max(R.sol.Z, [], 1);
    case {'hard-em', 'gibbs'}
        Zm = R.sol.Z;
end
visualize_results(K, X, U, Cx, Zm);


%% Sub functions


function visualize_results(K, X, U, Cx, Zm)

Gs = gaussd('m', U, Cx);

gm = intgroup(K, Zm);

hfig = figure;
set(hfig, 'Name', 'GMM Demo');

set(hfig, 'Position', [0 0 1024, 512]);
movegui(hfig, 'center');

subplot('Position', [0.05, 0.1, 0.43, 0.8]);
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 5);
hold on;
gaussd_ellipse(Gs, 2, 500, 'r-', 'LineWidth', 2);
axis equal;

subplot('Position', [0.52, 0.1, 0.43, 0.8]);

colors = {'r', 'g', 'b', 'm', 'c', 'k'};
for k = 1 : K   
    cr = colors{mod(k-1, numel(colors)) + 1};
    hold on;
    plot(X(1, gm{k}), X(2, gm{k}), [cr '.'], 'MarkerSize', 8);
end

axis equal;

