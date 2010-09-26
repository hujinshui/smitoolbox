function gmm_demo()
% A demo of using the facilities in smitoolbox to write gmm

% data configuration

n = [1000 2000 1000];
d = 2;
K = length(n);
alpha = 500;

mu0 = randn(d, K) * 3;
sigma = pdmat_udiag(d, 1);

g0_mp = gaussd_mp(mu0, sigma);
gpri = gaussd_mp(zeros(d, 1), pdmat_udiag(d, 1e8));

% produce data

sampler = g0_mp.get_sampler();
X = sampler.draw(1:K, n);

% show data

figure;
plot(X(1,:), X(2,:), '.');
draw_gauss(g0_mp, 'r-');

title('Ground truth');


% Initialize Qmap

Q0 = l2mat(K, repnum(n));
Q0 = Q0 + rand(K, sum(n)) * 0.2;
Q0 = bsxfun(@times, Q0, 1 ./ sum(Q0, 1));

% prepare estimation algorithm

gmm_est = fmm_est();

dpri = dirichletd(alpha * ones(K, 1));
qlabeler = multlabeler_std(dpri);
modelest = gaussd_map(gpri, sigma);

gmm_est.qlabeler = qlabeler;
gmm_est.model_est = modelest;
gmm_est.maxiter = 500;

% run estimation

problem = gmm_est.accept(X, Q0);
gmm_mon = fmm_est_mon(problem); %#ok<NASGU>  % create a monitor of the procedure

problem.solve();

% show result

figure;
plot(X(1,:), X(2,:), '.');
g_r = gaussd_mp(problem.models, sigma);
draw_gauss(g_r, 'r-');
title('Estimated models');


%% subfunctions

function draw_gauss(gs, varargin)

T = chol(gs.sigma);

t = linspace(0, 2*pi, 300);
circ = [cos(t); sin(t)];

for k = 1 : gs.nmodels
    
    if T.num == 1
        Y = bsxfun(@plus, T.take(1) * circ, gs.mu(:,k));
    else
        Y = bsxfun(@plus, T.take(k) * circ, gs.mu(:,k));
    end
    
    hold on;
    plot(Y(1,:), Y(2,:), varargin{:});
end







