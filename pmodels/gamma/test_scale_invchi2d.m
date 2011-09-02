function test_scale_invchi2d()
% Unit testing for scale_invchi2d class
%
%   test_scale_invchi2d;
%

%% main skeleton

do_test_on(1, false);
do_test_on(1, true);
do_test_on(3, false);
do_test_on(3, true);


%% core testing function

function do_test_on(d, is_uniform)

fprintf('Test on d = %d, is_uniform = %d ...\n', d, is_uniform);

nu = rand() + 8;

if is_uniform
    sigma2 = rand() + 2;
else
    sigma2 = rand(d, 1) + 2;
end

% test construction

if is_uniform
    D = scale_invchi2d(nu, sigma2, d);
    Sig2 = sigma2(ones(d, 1), 1);
else
    D = scale_invchi2d(nu, sigma2);
    Sig2 = sigma2;
end

% basic verification

assert(D.num == 1);
assert(D.dim == d);
assert(isequal(D.nu, nu));
assert(isequal(D.sigma2, sigma2));

alpha0 = nu / 2;
beta0 = nu * sigma2 / 2;

assert(norm(D.alpha - alpha0) < 2e-15);
assert(norm(D.beta - beta0) < 2e-15);

% verification of statistics

mean0 = (nu * sigma2) / (nu - 2);
var0 = (2 * nu^2 .* (sigma2.^2)) / ((nu - 2)^2 * (nu - 4));
mode0 = (nu * sigma2) / (nu + 2);
if is_uniform
    ent0 = d * (invgamma_entropy(D.alpha) + log(D.beta));
else
    ent0 = d * invgamma_entropy(D.alpha) + sum(log(D.beta));
end

assert(isequal(size(mean(D)), [d 1]));
assert(isequal(size(var(D)), [d 1]));
assert(isequal(size(mode(D)), [d 1]));
assert(isscalar(entropy(D)));

devcheck('mean calc', mean(D), mean0, 1e-13);
devcheck('var calc', var(D), var0, 1e-13);
devcheck('mode calc', mode(D), mode0, 1e-13);
devcheck('entropy', entropy(D), ent0, 1e-13);

% verification of logpdf and pdf

n = 100;
X = rand(d, n) + 1;

L0 = my_calc_logpdf(nu, Sig2, X);
L = D.logpdf(X);
assert(isequal(size(L), [1, n]));
devcheck('logpdf', L, L0, 1e-12);

P = D.pdf(X);
assert(isequal(size(P), [1, n]));
devcheck('pdf', P, exp(L), 1e-13);

% verification of sampling

ns = 5e6;
X = D.sample(ns);

devcheck('sample mean', vecmean(X), mean(D), 2e-2);
devcheck('sample var', vecvar(X), var(D), 5e-2);


%% Auxiliary functions

function L = my_calc_logpdf(nu, sigma2, X)

d = size(sigma2, 1);
L = zeros(d, size(X,2));
for i = 1 : d
    x = X(i, :);
    c = (nu / 2) * log(sigma2(i) * nu / 2) - gammaln(nu / 2);
    t1 = (nu * sigma2(i)) ./ (2 * x);
    t2 = (1 + nu / 2) * log(x);
    L(i, :) = c - t1 - t2;
end
L = sum(L, 1);


function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_scale_invchi2d:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end



