function test_gaussd(pdmatclassname, cfgs, est)
% A unit-test function to test the correctness of Gaussian models
%
%   test_gaussd(pdmatclassname, cfgs);
%   test_gaussd(pdmatclassname, cfgs, est);
%       Tests the correctness of the implementation of Gaussian models.
%
%       pdmatclassname is the name of the class for representing the
%       covariance matrix. 
%   
%       cfgs is an n x 3 array, where each row gives the configuration
%       for a batch of testing in form of [d, m, is_shared]
%
%       est is an MLE estimator. If specified, MLE estimation is also
%       tested.
%

%% main skeleton

assert(ischar(pdmatclassname) && ismember('pdmat', superclasses(pdmatclassname)), ...
    'test_gaussd:invalidarg', ...
    'pdmatclassname should be the name of a derived class of pdmat.');

if nargin < 3
    est = [];
end

if ~isempty(est)
    assert(isa(est, 'gmle_base'));
end

ncfgs = size(cfgs, 1);
for i = 1 : ncfgs
    test_on_config(pdmatclassname, cfgs(i,1), cfgs(i,2), cfgs(i,3), est);
end
        

%% main test function

function test_on_config(gmname, d, m, is_shared, est)

% generate random models

if is_shared
    m2 = 1;
else
    m2 = m;
end

mu = randn(d, m);
sigma = feval([gmname '.random'], d, m2);

g_cp = gaussd_cp.from_mean_and_cov(mu, sigma);
g_mp = gaussd_mp(mu, sigma);

% basic verification

assert(g_cp.dim == d);
assert(g_cp.nmodels == m);

assert(g_mp.dim == d);
assert(g_mp.nmodels == m);

fprintf('Test on pdmat-class: %s [d = %d, m = %d, s = %d]\n', gmname, d, m, is_shared);

% Test parameter consistency

test_param_consistency(g_cp, g_mp, d, m, is_shared);

% Test conversion

test_cp_mp_conv(g_cp, g_mp);

% Test evaluation

test_eval(g_cp);

% Test sampling and MLE estimation

if ~isempty(est)    
    test_sampler_and_mle(g_mp, est);
end



%% Core testing functions

function test_param_consistency(g_cp, g_mp, d, m, is_shared)

assert(isequal(size(g_cp.theta1), [d m]));
assert(g_cp.theta2.dim == d);
if is_shared
    assert(g_cp.theta2.num == 1);
else
    assert(g_cp.theta2.num == m);
end
assert(isequal(size(g_cp.theta3), [1 m]));

assert(isequal(size(g_mp.mu),[d m]));
assert(g_mp.sigma.dim == d);
if is_shared
    assert(g_mp.sigma.num == 1);
else
    assert(g_mp.sigma.num == m);
end

inv_theta2 = inv(g_cp.theta2);
if is_shared
    theta3_r = inv_theta2.quadterm(g_cp.theta1);
else
    theta3_r = sum(g_cp.theta1 .* cmult(inv_theta2, g_cp.theta1), 1); 
end
devcheck('theta3-consist', g_cp.theta3, theta3_r, 1e-10);

logpar_r = 0.5 * (d * log(2*pi) - g_cp.theta2.logdet + theta3_r);
devcheck('logpar', g_cp.logpar, logpar_r, 1e-10);


function test_cp_mp_conv(g_cp, g_mp)

g_mp2 = to_mp(g_cp);
g_cp2 = to_cp(g_mp);

devcheck('compare-mp-mp2-mu', g_mp.mu, g_mp2.mu, 1e-10);
devcheck('compare-mp-mp2-sigma', g_mp.sigma.fullform, g_mp2.sigma.fullform, 1e-10);

devcheck('compare-cp-cp2-theta1', g_cp.theta1, g_cp2.theta1, 1e-10);
devcheck('compare-cp-cp2-theta2', g_cp.theta2.fullform, g_cp2.theta2.fullform, 1e-10);
devcheck('compare-cp-cp3-theta3', g_cp.theta3, g_cp2.theta3, 1e-10);


function test_eval(g)

d = g.dim;
m = g.nmodels;

T1 = g.theta1;
T2 = g.theta2.fullform;
mu = zeros(d, m);
ld = zeros(1, m);

if size(T2, 3) == 1 && m > 1
    T2 = repmat(T2, [1, 1, m]);
end

for i = 1 : m   
    mu(:,i) = T2(:,:,i) \ T1(:,i);
    ld(i) = calc_logdet(T2(:,:,i));
end

n0 = 100;
X0 = randn(d, n0);

Dsq0 = comp_sq_mahdist(mu, T2, X0);
Dsq1 = mahdist_sq(g, X0);
Dsq2 = zeros(m, n0);
for i = 1 : m
    Dsq2(i, :) = mahdist_sq(g, X0, i);
end

D0 = sqrt(max(Dsq0, 0));
D1 = mahdist(g, X0);
D2 = zeros(m, n0);
for i = 1 : m
    D2(i, :) = mahdist(g, X0, i);
end

assert(isequal(size(Dsq1), size(D1), [m, n0]));
devcheck('mahdist_sq_1', Dsq0, Dsq1, 1e-8 * max(Dsq0(:)));
devcheck('mahdist_sq_2', Dsq0, Dsq2, 1e-8 * max(Dsq0(:)));
devcheck('mahdist_1', D0, D1, 1e-8 * max(D0(:)));
devcheck('mahdist_2', D0, D2, 1e-8 * max(D0(:)));

LP0 = bsxfun(@plus, -0.5 * Dsq0, 0.5 * (ld.' - d * log(2*pi)));
LP1 = logprob(g, X0);
LP2 = zeros(m, n0);
for i = 1 : m
    LP2(i, :) = logprob(g, X0, i);
end

assert(isequal(size(LP1), [m n0]));
devcheck('logpdf_1', LP0, LP1, 1e-8 * max(abs(LP0(:))));
devcheck('logpdf_2', LP0, LP2, 1e-8 * max(abs(LP0(:))));

if d == 1
    pint = arrayfun(@(i) prob_integrate_1d(g, i, mu(:,i), T2(:,:,i)), 1:m);
    devcheck('prob_integrate_1d', pint, ones(1, m), 1e-7); 
elseif d == 2
    pint = arrayfun(@(i) prob_integrate_2d(g, i, mu(:,i), T2(:,:,i)), 1:m);
    devcheck('prob_integrate_2d', pint, ones(1, m), 1e-5);
end


function test_sampler_and_mle(g, est)

assert(isa(g, 'gaussd_mp'));
n = 100000;

d = g.dim;
m = g.nmodels;

spl = g.get_sampler();
X = spl.draw([], n);

assert(isequal(size(X), [d, m * n]));

ws = cell(1, m);
for i = 1 : m
    ws{i} = ones(1, n);
end
W = blkdiag(ws{:});

est.tiec = g.sigma.num == 1;
est.to_cp = false;

g1 = est.estimate(X, W);
assert(isa(g1, 'gaussd_mp'));
assert(g1.dim == g.dim && g1.nmodels == g.nmodels);
assert(isequal(size(g.mu), size(g1.mu)));
assert(isa(g1.sigma, class(g.sigma)));
assert(g1.sigma.dim == g.sigma.dim && g1.sigma.num == g.sigma.num); 

devcheck('est_mu', g.mu, g1.mu, 1e-2);
devcheck('est_sigma', g.sigma.fullform, g1.sigma.fullform, 1e-2);

if m == 1
    g2 = est.estimate(X);
    
    assert(isequal(size(g.mu), size(g2.mu)));
    assert(isequal(size(g.sigma), size(g2.sigma)));
    
    devcheck('est_mu_noweight', g1.mu, g2.mu, 1e-8);
    devcheck('est_sigma_noweight', g1.sigma.fullform, g2.sigma.fullform, 1e-8);
end



%% Auxiliary function

function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_gaussd:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end


function v = calc_logdet(C)

L = chol(C);
v = 2 * sum(log(diag(L)));


function D = comp_sq_mahdist(mu, icov, X)
% compute squared Mahalanobis distance to Gaussian centers

m = size(mu, 2);
n = size(X, 2);

D = zeros(m, n);
for i = 1 : m    
    A = icov(:,:,i);
    v = mu(:,i);
    D(i, :) = pwmahdist(v, X, A, 'square');
end

function v = prob_integrate_1d(g, i, mu, icov)

assert(g.dim == 1);
s = sqrt(1 / icov);

n = 2000;
x = linspace(mu - 10 * s, mu + 10 * s, n);
p = g.prob(x, i);

v = trapz(x, p);

function v = prob_integrate_2d(g, i, mu, icov)

assert(g.dim == 2);
C = inv(icov);

[V, D] = eig(C);
v1 = V(:,1);
v2 = V(:,2);
s1 = sqrt(D(1,1));
s2 = sqrt(D(2,2));

[X1, X2] = meshgrid(-8:0.1:8, -8:0.1:8);

X = bsxfun(@times, X1(:)' * s1, v1) + bsxfun(@times, X2(:)' * s2, v2);
X = bsxfun(@plus, X, mu);

p = g.prob(X, i);

v = sum(p) * (s1 * s2) * 0.01;



