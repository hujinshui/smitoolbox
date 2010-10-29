function rgmm_demo()
% A demo of using the facilities in smitoolbox to write robust gmm

% data configuration

n = [1000 2000 1000];
n_out = 1000;
d = 2;
K = length(n);
alpha = 500;
b_alpha = 400;
b_beta = 100;
lpout = -5;

mu0 = randn(d, K) * 3;
sigma = pdmat_udiag(d, 1);

g0_mp = gaussd_mp(mu0, sigma);
gpri = gaussd_mp(zeros(d, 1), pdmat_udiag(d, 1e8));

% produce data

sampler = g0_mp.get_sampler();
X1 = sampler.draw(1:K, n);
X0 = 20 * (rand(d, n_out) - 0.5);

X = [X1 X0];

% show data

figure;
plot(X(1,:), X(2,:), '.');
draw_gauss(g0_mp, 'r-');

title('Ground truth');


% Initialize Qmap and Gmap

Q0 = l2mat(K, repnum(n));
Q0 = [Q0 + rand(K, sum(n)) * 0.2, rand(K, n_out)];
Q0 = bsxfun(@times, Q0, 1 ./ sum(Q0, 1));

G0 = 0.9;

% prepare estimation algorithm

rgmm_est = rfmm_est();

dpri = dirichletd(alpha * ones(K, 1));
bpri = betad(b_alpha, b_beta);

qlabeler = multlabeler_std(dpri);
glabeler = binlabeler_std(bpri);
modelest = gaussd_map(gpri, sigma);

rgmm_est.qlabeler = qlabeler;
rgmm_est.glabeler = glabeler;
rgmm_est.model_est = modelest;
rgmm_est.lp_outlier = lpout;

rgmm_est.maxiter = 500;

% run estimation

problem = rgmm_est.accept(X, Q0, G0);
rgmm_mon = rfmm_est_mon(problem); %#ok<NASGU>  % create a monitor of the procedure

problem.solve();

% show result

figure;

g = problem.gposterior;
scatter(X(1,:), X(2,:), 1 + 5 * g);
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







