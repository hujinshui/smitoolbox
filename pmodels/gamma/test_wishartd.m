function test_wishartd()
% Unit testing of wishartd class and wishart_sample function
%
%   test_wishartd();
%

%% main

ds = [1 2 3 5];
tys = {'s', 'd', 'f'};

for d = ds
    for it = 1 : numel(tys)
        do_test_on(d, tys{it});
    end
end

%%  core testing function

function do_test_on(d, ty)

fprintf('test on dim = %d, ty = %s\n', d, ty);

% construction

S = randpdm(ty, d, 1);
deg = 8.5; 
Sm = pdmat_fullform(S);

D0 = wishartd(S, deg);
D1 = wishartd(S, deg, 'pre');

assert(D0.num == 1);
assert(D0.dim == d);
assert(D0.deg == deg);
assert(isequal(D0.S, S));
assert(isempty(D0.invS));
assert(isempty(D0.lpconst));

assert(D1.num == 1);
assert(D1.dim == d);
assert(D1.deg == deg);
assert(isequal(D1.S, S));
assert(~isempty(D1.invS));
assert(~isempty(D1.lpconst));

% verify pre-computation

q = deg / 2;
lpc = q*d*log(2) + q * pdmat_lndet(S) + mvgammaln(d, q); 

devcheck('lpconst', D1.lpconst, lpc, 1e-13);
assert(isequal(pdmat_inv(S), D1.invS));

% verify statistics

Mean0 = deg * Sm;
Mode0 = (deg - d - 1) * Sm;

devcheck('mean (D0)', pdmat_fullform(mean(D0)), Mean0, 1e-13);
devcheck('mean (D1)', pdmat_fullform(mean(D1)), Mean0, 1e-13);
devcheck('mode (D0)', pdmat_fullform(mode(D0)), Mode0, 1e-13);
devcheck('mode (D1)', pdmat_fullform(mode(D1)), Mode0, 1e-13);

% verify logpdf evaluation

n = 100;
W_s = randpdm('s', d, n);
W_d = randpdm('d', d, n);
W_f = randpdm('f', d, n);

L1_s0 = D0.logpdf(W_s);
L1_d0 = D0.logpdf(W_d);
L1_f0 = D0.logpdf(W_f);

L1_s1 = D1.logpdf(W_s);
L1_d1 = D1.logpdf(W_d);
L1_f1 = D1.logpdf(W_f);

assert(isequal(size(L1_s0), [1, n]));
assert(isequal(size(L1_d0), [1, n]));
assert(isequal(size(L1_f0), [1, n]));

assert(isequal(L1_s0, L1_s1));
assert(isequal(L1_d0, L1_d1));
assert(isequal(L1_f0, L1_f1));

L0_s = my_logpdf(Sm, deg, W_s, lpc);
L0_d = my_logpdf(Sm, deg, W_d, lpc);
L0_f = my_logpdf(Sm, deg, W_f, lpc);

devcheck('loglik (s)', L1_s0, L0_s, 1e-12);
devcheck('loglik (d)', L1_d0, L0_d, 1e-12);
devcheck('loglik (f)', L1_f0, L0_f, 1e-12);

% verify standard sample

N0 = 50;
N = 2e5;

Y = wishart_sample(d, deg);
assert(isequal(size(Y), [d d]));
assert(isequal(Y, Y'));
assert(all(eig(Y)) > 0);

W = D0.sample();
assert(is_pdmat(W) && W.d == d && W.n == 1);
assert(W.ty == 'f');
assert(isequal(W.v, W.v'));
assert(all(eig(W.v)) > 0);

Ys = wishart_sample(d, deg, N);
assert(isequal(size(Ys), [d d N]));
for i = 1 : N0
    Yi = Ys(:,:,i);
    assert(isequal(Yi, Yi'));
    assert(all(eig(Yi)) > 0);
end

Ws = D0.sample(N);
assert(is_pdmat(Ws) && Ws.d == d && Ws.n == N);
for i = 1 : N0
    Wi = pdmat_fullform(Ws, i);
    assert(isequal(Wi, Wi'));
    assert(all(eig(Wi)) > 0);
end

Ymean0 = deg * eye(d);
Wmean0 = deg * Sm;

Id = eye(d);
Yvar0 = zeros(d, d);
Wvar0 = zeros(d, d);
for i = 1 : d
    for j = 1 : d
        Yvar0(i, j) = deg * (Id(i,j)^2 + Id(i,i) * Id(j,j));
        Wvar0(i, j) = deg * (Sm(i,j)^2 + Sm(i,i) * Sm(j,j));
    end
end

Ymean = mean(Ys, 3);
Yvar = reshape(vecvar(reshape(Ys, d*d, N)), d, d);
Wmean = mean(Ws.v, 3);
Wvar = reshape(vecvar(reshape(Ws.v, d*d, N)), d, d);

devcheck('sample mean (std)', Ymean, Ymean0, 5e-2);
devcheck('sample var (std)', Yvar, Yvar0, 0.3);
devcheck('sample mean', Wmean, Wmean0, 5e-2);
devcheck('sample mean', Wvar, Wvar0, 0.3);


%% Auxiliary functions

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


function L = my_logpdf(Sm, m, W, lpc)

d = W.d;
n = W.n;
L = zeros(1, n);

for i = 1 : n
    Wi = pdmat_fullform(W, i);
    
    u1 = (m - d - 1) / 2 * lndet(Wi);
    u2 = (-1/2) * trace(Sm \ Wi);
    L(i) = u1 + u2 - lpc;
end


function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_wishartd:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end


