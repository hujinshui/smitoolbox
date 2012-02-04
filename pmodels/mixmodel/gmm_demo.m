function R = gmm_demo(cf, op)
% A script to demo the inference over Gaussian mixture model
%
%   R = gmm_demo(cf);
%   R = gmm_demo(cf);
%
%       Runs a Gaussian mixture model demo using the specified covariance
%       matrix form, which can be either 's' (scale form), 
%       'd' (diagonal form), or 'f' (full matrix form).
%
%       The output is the struct that captures inference result.
%
%   R = gmm_demo(cf, 'tied-cov');
%   R = gmm_demo(cf, 'tied-cov');
%
%       Runs the demo under the setting that the covariance matrix
%       is shared across all mixture components.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 14, 2010
%       - Modified by Dahua Lin, on Aug 31, 2011
%       - Modified by Dahua Lin, on Sep 28, 2011
%       - Modified by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~(ischar(cf) && isscalar(cf) && (cf == 's' || cf == 'd' || cf == 'f'))
    error('gmm_demo:invalidarg', 'Invalid cf.');
end

if nargin < 2
    c_tied = false;
    mitrs = 200;
else
    if ~(ischar(op) && strcmpi(op, 'tied-cov'))
        error('gmm_demo:invalidarg', 'The second argument is invalid.');
    end
    c_tied = true;
    mitrs = 1600;
end


%% generate data

d = 2;
K = 3;
n = 1000;   % #samples / class

pri_sigma = 20;
gpri = gaussd('m', 0, pdmat('s', 2, pri_sigma));

U0 = gaussd_sample(gpri, K);

Xs = cell(1, K);

if c_tied
    Cx = rand_pdmat(cf, d, 1, [0.5 1.5]);    
    for k = 1 : K
        gx = gaussd('m', U0(:,k), Cx);
        Xs{k} = gaussd_sample(gx, n);
    end
else
    for k = 1 : K
        Cx = rand_pdmat(cf, d, 1, [0.5, 1.5]);
        gx = gaussd('m', U0(:,k), Cx);
        Xs{k} = gaussd_sample(gx, n);
    end
end
X = [Xs{:}];

%% Run estimation

w = [];
R = gmm_fit(X, w, K, ...
    'cov_form', cf, 'pricount', 0, 'tied_cov', c_tied, ...
    'maxiters', mitrs, 'tol', 1e-6, 'display', 'eval');

%% Visualize

[~, Zm] = max(R.Q, [], 1);
visualize_results(K, X, R.G, Zm);



%% Sub functions

function visualize_results(K, X, Gs, Zm)

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



