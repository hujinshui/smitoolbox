function test_gaussd()
% A unit-test function to test the correctness of gaussd implementation
%
%   test_gaussd;
%       Tests the correctness of implementation of gaussd with a 
%       specified covariance matrix class.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jun 19, 2010
%       - Modified by Dahua Lin, on Aug 17, 2011
%       - Modified by Dahua Lin, on Aug 25, 2011
%

%% main skeleton

cfs = {'s', 'd', 'f'};
ds = [1, 2, 5];
ms = [1, 3];

for icf = 1 : numel(cfs)
    cf = cfs{icf};
    for d = ds
        for m = ms
            test_on_config(cf, d, m, 0, 0);
            test_on_config(cf, d, m, 0, 1);
            test_on_config(cf, d, m, 1, 0);
        end
    end
end
        

%% main test function

function test_on_config(cf, d, m, is_zeromean, is_shared)

if is_zeromean && m > 1
    return;
end

fprintf('Test gaussd for with cf = %c, dim = %d, num = %d, zmean = %d, sharedcov = %d\n', ...
    cf, d, m, is_zeromean, is_shared);

% generate random models

if is_shared
    m2 = 1;
else
    m2 = m;
end

if is_zeromean
    mu = zeros(d, m);
else
    mu = randn(d, m);
end
    
[sigma, vars] = randpdm(cf, d, m2);
if is_shared
    vars = repmat(vars, 1, m);
end

g_mp0 = gaussd.from_mp(mu, sigma);
g_mp1 = gaussd.from_mp(mu, sigma, 'ip');

assert(g_mp0.zmean == is_zeromean);
assert(g_mp0.shared_cov == (m2 == 1));
assert(g_mp1.zmean == is_zeromean);
assert(g_mp1.shared_cov == (m2 == 1));

h = g_mp1.h;
J = g_mp1.J;
c0 = g_mp1.c0;

g_cp0 = gaussd.from_ip(h, J, c0);
g_cp1 = gaussd.from_ip(h, J, c0, 'mp');

g_cp0_a = gaussd.from_ip(h, J);
g_cp1_a = gaussd.from_ip(h, J, [], 'mp');

assert(g_cp0.zmean == is_zeromean);
assert(g_cp0.shared_cov == (m2 == 1));
assert(g_cp1.zmean == is_zeromean);
assert(g_cp1.shared_cov == (m2 == 1));

% basic verification

basic_verify(g_mp0, d, m, true, false);
basic_verify(g_mp1, d, m, true, true);
basic_verify(g_cp0, d, m, false, true);
basic_verify(g_cp1, d, m, true, true);
basic_verify(g_cp0_a, d, m, false, true);
basic_verify(g_cp1_a, d, m, true, true);

verify_mp(g_mp0, mu, sigma);
verify_mp(g_mp1, mu, sigma);
verify_cp(g_cp0, h, J, c0);
verify_cp(g_cp1, h, J, c0);
verify_cp(g_cp0_a, h, J);
verify_cp(g_cp1_a, h, J);

devcheck('c0 consistency', g_cp0_a.c0, c0, 1e-12);
devcheck('c0 consistency', g_cp1_a.c0, c0, 1e-12);

devcheck('mu consistency', g_cp1.mu, mu, 1e-12);
devcheck('sigma consistency', g_cp1.C.v, sigma.v, 1e-12);

% verification of statistics

assert(isequal(mean(g_mp0), mu));
assert(isequal(var(g_mp0), vars));
if m == 1 || is_shared
    assert(isequal(cov(g_mp0), pdmat_fullform(sigma)));
end
for i = 1 : m
    if is_shared
        assert(isequal(cov(g_mp0, i), pdmat_fullform(sigma)));
    else
        assert(isequal(cov(g_mp0, i), pdmat_fullform(sigma, i)));
    end
end

ents0 = zeros(1, m);
for i = 1 : m
    if is_shared
        sig_i = sigma;
    else
        sig_i = pdmat_pick(sigma, i);
    end
    ents0(i) = (1/2) * ( d * log( (2*pi*exp(1)) ) + pdmat_lndet(sig_i) );
end
ents1 = entropy(g_mp0);
ents2 = entropy(g_cp0);

ents1a = zeros(1, m);
ents2a = zeros(1, m);
for i = 1 : m
    ents1a(i) = entropy(g_mp0, i);
    ents2a(i) = entropy(g_cp0, i);
end

assert(isequal(size(ents1), [1, m]));
assert(isequal(size(ents2), [1, m]));

devcheck('entropy calc (C)', ents0, ents1, 1e-12);
devcheck('entropy calc (J)', ents0, ents2, 1e-12);
devcheck('entropy calc (C-i)', ents0, ents1a, 1e-12);
devcheck('entropy calc (J-i)', ents0, ents2a, 1e-12);

% Test evaluation
test_eval(g_cp1);

% Test sampling
test_sampling(g_mp1);


%% Evaluation testing

function test_eval(g)

d = g.dim;
m = g.num;

mu = g.mu;
if isequal(mu, 0)
    mu = zeros(d, m);
end

if m == 1
    C2 = pdmat_fullform(g.C);
else
    C2 = zeros(d, d, m);
    if g.shared_cov
        sC2 = pdmat_fullform(g.C);
        C2 = repmat(sC2, [1, 1, m]);
    else
        for k = 1 : m
            C2(:,:,k) = pdmat_fullform(g.C, k);
        end
    end
end

ld = zeros(1, m);
Q2 = zeros(d, d, m);
for k = 1 : m
    ld(k) = lndet(C2(:,:,k));
    Q2(:,:,k) = inv(C2(:,:,k));
end

n0 = 100;
X0 = randn(d, n0);

Dsq0 = comp_sq_mahdist(mu, Q2, X0);

Dsq1 = sqmahdist(g, X0);
Dsq2 = zeros(m, n0);
for i = 1 : m
    Dsq2(i, :) = sqmahdist(g, X0, i);
end

assert(isequal(size(Dsq1), [m n0]));
devcheck('sqmahdist_1', Dsq0, Dsq1, 1e-8 * max(abs(Dsq0(:))));
devcheck('sqmahdist_2', Dsq0, Dsq2, 1e-8 * max(abs(Dsq0(:))));

LP0 = bsxfun(@plus, -0.5 * Dsq0, 0.5 * (- ld.' - d * log(2*pi)));
LP1 = logpdf(g, X0);
LP2 = zeros(m, n0);
for i = 1 : m
    LP2(i, :) = logpdf(g, X0, i);
end

assert(isequal(size(LP1), [m n0]));
devcheck('logpdf_1', LP0, LP1, 1e-8 * max(abs(LP0(:))));
devcheck('logpdf_2', LP0, LP2, 1e-8 * max(abs(LP0(:))));

if d == 1
    pint = arrayfun(@(i) prob_integrate_1d(g, i, mu(:,i), Q2(:,:,i)), 1:m);
    devcheck('prob_integrate_1d', pint, ones(1, m), 1e-7); 
elseif d == 2
    pint = arrayfun(@(i) prob_integrate_2d(g, i, mu(:,i), Q2(:,:,i)), 1:m);
    devcheck('prob_integrate_2d', pint, ones(1, m), 1e-5);
end


function test_sampling(g)

d = g.dim;
m = g.num;
n = 100000;
ns = n * ones(1, m);

mu = g.mu;
if m == 1
    C = pdmat_fullform(g.C);
else
    C = zeros(d, d, m);
    if g.shared_cov
        sC = pdmat_fullform(g.C);
        C = repmat(sC, [1, 1, m]);
    else
        for k = 1 : m
            C(:,:,k) = pdmat_fullform(g.C, k);
        end
    end
end

if m == 1
    X = g.sample(ns);
else
    X = g.sample(ns, 1:m);
end
    
assert(size(X, 2) == sum(ns));

mu_s = zeros(d, m);
C_s = zeros(d, d, m);

for k = 1 : m    
    Xk = X(:, n * (k-1) + (1:n));
    
    c_mu = sum(Xk, 2) * (1/n);
    Zk = bsxfun(@minus, Xk, c_mu);
    c_sig = (Zk * Zk') * (1/n);
    
    mu_s(:,k) = c_mu;
    C_s(:,:,k) = c_sig;
end

devcheck('sample_mean', mu, mu_s, 2e-2);
devcheck('sample_cov', C, C_s, 5e-2);


%% Auxiliary function

function [C, vs] = randpdm(cf, d, m)

switch cf
    case 's'
        v = 0.5 + rand(1, m);
        C = pdmat(cf, d, v);
        vs = repmat(v, d, 1);
    case 'd'
        v = 0.5 + rand(d, m);
        C = pdmat(cf, d, v);
        vs = v;
    case 'f'
        C = zeros(d, d, m);
        for i = 1 : m
            R = orth(randn(d));
            dv = 0.5 + rand(d, 1);
            C(:,:,i) = R' * diag(dv) * R;
        end
        C = pdmat(cf, d, C);
        vs = zeros(d, m);
        for i = 1 : m
            vs(:, i) = diag(C.v(:,:,i));
        end
end


function basic_verify(g, d, m, has_mp, has_ip)

assert(g.dim == d && g.num == m && g.has_mp == has_mp && g.has_ip == has_ip);


function verify_mp(g, mu, sigma)

if size(mu,2) == 1 && all(mu == 0)
    assert(isequal(g.mu, 0) && g.zmean)
else
    assert(isequal(g.mu, mu) && ~g.zmean);
end
assert(isequal(g.C, sigma));


function verify_cp(g, h, J, c0)

assert(isequal(g.h, h) && isequal(g.J, J));

if nargin >= 4
    assert(isequal(g.c0, c0));
end


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
    D(i, :) = pwsqmahdist(v, X, A);
end

function v = prob_integrate_1d(g, i, mu, icov)

assert(g.dim == 1);
s = sqrt(1 / icov);

n = 2000;
x = linspace(mu - 10 * s, mu + 10 * s, n);
p = g.pdf(x, i);

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

p = g.pdf(X, i);

v = sum(p) * (s1 * s2) * 0.01;



