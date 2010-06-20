function test_gaussd(symatclassname, ds, ms)
% A unit-test function to test the correctness of gaussd implementation
%
%   test_gaussd(symatclassname, ds, ms);
%       Tests the correctness of implementation of gaussd with a 
%       specified covariance matrix class.
%
%       In the input, symatclassname is the name of the class for
%       representing the covariance matrix. ds and ms are vectors
%       of dimensions and model numbers for testing.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jun 19, 2010
%

%% main skeleton

for d = ds
    for m = ms
        test_on_config(symatclassname, d, m, 0, 0);
        test_on_config(symatclassname, d, m, 0, 1);
        test_on_config(symatclassname, d, m, 1, 0);
    end
end
        

%% main test function

function test_on_config(symatclassname, d, m, is_zeromean, is_shared)

fprintf('Test gaussd for [%s] with dim = %d, num = %d, zmean = %d, sharedcov = %d\n', ...
    symatclassname, d, m, is_zeromean, is_shared);

% generate random models

if is_shared
    m2 = 1;
else
    m2 = m;
end

if is_zeromean
    mu = 0;
else
    mu = randn(d, m);
end
    
sigma = feval([symatclassname '.randpdm'], d, m2);

G0 = gaussd.from_mp(mu, sigma);
G = gaussd.from_mp(mu, sigma, 'cp');
g_cp = gaussd.from_cp(G.c1, G.c2);
g_cp2 = gaussd.from_cp(G.c1, G.c2, G.c0);

% basic verification

assert(G0.dim == d && G0.num == m && G0.use_mp && ~G0.use_cp);
assert(isequal(G0.mu, G.mu, mu));
assert(isequal(fullform(G0.sigma), fullform(G.sigma), fullform(sigma)));

assert(g_cp.dim == d && g_cp.num == m && g_cp.use_mp && g_cp.use_cp);
assert(g_cp2.dim == d && g_cp2.num == m && g_cp2.use_mp && g_cp2.use_cp);
assert(isequal(G.c1, g_cp.c1, g_cp2.c1));
assert(isequal(G.c2, g_cp.c2, g_cp2.c2));
assert(isequal(G.c0, g_cp2.c0));
devcheck('c0 consistency', G.c0, g_cp.c0, 1e-12);
assert(isequal(g_cp.mu, g_cp2.mu));
assert(isequal(fullform(g_cp.sigma), fullform(g_cp2.sigma)));

devcheck('mu consistency', G.mu, g_cp.mu, 1e-12);
devcheck('sigma consistency', fullform(G.sigma), fullform(g_cp.sigma), 1e-12);

% Test evaluation

test_eval(G);

%% Evaluation testing

function test_eval(g)

d = g.dim;
m = g.num;

mu = g.mu;
if isequal(mu, 0)
    mu = zeros(d, m);
end
ld = lndet(g.sigma);
T2 = fullform(g.c2);
if ndims(T2) == 2 && m > 1
    T2 = repmat(T2, [1, 1, m]);
end

n0 = 100;
X0 = randn(d, n0);

Dsq0 = comp_sq_mahdist(mu, T2, X0);
LP0 = bsxfun(@plus, -0.5 * Dsq0, 0.5 * (- ld.' - d * log(2*pi)));
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


%% Auxiliary function

function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_gaussd:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end

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
p = exp(g.logprob(x, i));

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

p = exp(g.logprob(X, i));

v = sum(p) * (s1 * s2) * 0.01;



