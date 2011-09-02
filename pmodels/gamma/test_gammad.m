function test_gammad()
% Unit testing for Gamma distribution
%
%   test_gammad;
%

% Create by Dahua Lin, on Sep 1, 2011
%

%% main skeleton

ds = [1 2 5];
ms = [1 3];

for d = ds
    for m = ms        
        do_test_on(d, m, 0, 0);
        do_test_on(d, m, 0, 1);
        do_test_on(d, m, 1, 0);
        do_test_on(d, m, 1, 1);
    end
end


%% core test function

function do_test_on(d, m, uniform_shape, share_scale)

fprintf('test on d = %d, m = %d, uniform_shape = %d, share_scale = %d ...\n', ...
    d, m, uniform_shape, share_scale);

% construction

if uniform_shape
    alpha = rand(1, m) + 1.2;
    A = repmat(alpha, [d, 1]);
else
    alpha = rand(d, m) + 1.2;
    A = alpha;
end

if share_scale
    beta = rand() + 0.5;
    B = repmat(beta, 1, m);
else
    beta = rand(1, m) + 0.5;
    B = beta;
end

if ~uniform_shape
    g =  gammad(alpha, beta);
    gp = gammad(alpha, beta, [], 'pre');
else
    g =  gammad(alpha, beta, d);
    gp = gammad(alpha, beta, d, 'pre');
end

% basic verification

assert(g.dim == d);
assert(g.num == m);
assert(isequal(g.alpha, alpha));
assert(isequal(g.beta, beta));
assert(isempty(g.lpconst));

assert(gp.dim == d);
assert(gp.num == m);
assert(isequal(gp.alpha, alpha));
assert(isequal(gp.beta, beta));
assert(isequal(size(gp.lpconst), [1, m]));

% verify statistics computation

[mean0, var0, mode0, ent0] = my_calc_stats(A, B);

mean1 = mean(g);
var1 = var(g);
mode1 = mode(g);
ent1 = entropy(g);

assert(isequal(size(mean1), [d, m]));
assert(isequal(size(var1), [d, m]));
assert(isequal(size(mode1), [d, m]));
assert(isequal(size(ent1), [1, m]));

devcheck('mean calc', mean0, mean1, 1e-14);
devcheck('var calc',  var0, var1, 1e-14);
devcheck('mode calc',  mode0, mode1, 1e-14);
devcheck('entropy calc', ent0, ent1, 1e-12);

% verify the computation of lpconst

lpc0 = my_calc_lpconst(A, B);
devcheck('lpconst calc', gp.lpconst, lpc0, 1e-12); 

% verify the computation of logpdf & pdf

N = 100;
X = rand(d, N);

L0 = my_calc_logpdf(A, B, X);
L1 = g.logpdf(X);
L1p = gp.logpdf(X);
assert(isequal(size(L1), [m, N]));
assert(isequal(L1, L1p));

L2 = zeros(m, N);
for k = 1 : m
    L2(k, :) = gp.logpdf(X, k);
end

devcheck('logpdf eval', L1, L0, 1e-10);
devcheck('logpdf eval (per-row)', L2, L1, 1e-13);

P1 = g.pdf(X);
assert(isequal(P1, exp(L1)));

% verify sampling

ns = 5e5;
if m == 1
    X1 = g.sample(ns);
    assert(isequal(size(X1), [d, ns]));
    
    devcheck('sample 1 - mean', vecmean(X1), mean1, 2e-2);
    devcheck('sample 1 - var',  vecvar(X1), var1, 5e-2);    
end

X2 = g.sample(ns(ones(1, m)), 1:m);
for k = 1 : m
    cX2 = X2(:, (k-1)*ns+1 : (k-1)*ns+ns);
        
    devcheck('sample 2 - mean', vecmean(cX2), mean1(:,k), 2e-2);
    devcheck('sample 2 - var',  vecvar(cX2), var1(:,k), 5e-2);   
end




%% Auxiliary functions

function v = my_calc_lpconst(A, B)

v = bsxfun(@times, A, log(B)) + gammaln(A);

if size(v, 1) > 1
    v = sum(v, 1);
end
m = size(A, 2);
assert(isequal(size(v), [1, m]));
v = -v;


function [M, V, Mo, E] = my_calc_stats(A, B)

M = bsxfun(@times, A, B);
V = bsxfun(@times, A, B.^2);
Mo = bsxfun(@times, A-1, B);

E1 = A + gammaln(A) + (1 - A) .* psi(A);
E = bsxfun(@plus, E1, log(B));
E = sum(E, 1);


function L = my_calc_logpdf(A, B, X)

[d, m] = size(A);
n = size(X, 2);

L = zeros(m, n);

for k = 1 : m
    
    b = B(k);
        
    Pk = zeros(d, n);
    for i = 1 : d
        a = A(i, k); 
        Pk(i, :) = gampdf(X(i, :), a, b);
    end
    
    L(k, :) = sum(log(Pk), 1);
end


function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_gammad:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end



