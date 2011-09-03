function test_invwishartd()
% Unit testing of invwishartd class 
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

Phi = randpdm(ty, d, 1);
deg = 8.5; 
Phim = pdmat_fullform(Phi);

D0 = invwishartd(Phi, deg);
D1 = invwishartd(Phi, deg, 'pre');

assert(D0.num == 1);
assert(D0.dim == d);
assert(D0.deg == deg);
assert(isequal(D0.Phi, Phi));
assert(isempty(D0.invPhi));
assert(isempty(D0.lpconst));

assert(D1.num == 1);
assert(D1.dim == d);
assert(D1.deg == deg);
assert(isequal(D1.Phi, Phi));
assert(~isempty(D1.invPhi));
assert(~isempty(D1.lpconst));

% verify pre-computation

q = deg / 2;
lpc = q*d*log(2) - q * pdmat_lndet(Phi) + mvgammaln(d, q); 

devcheck('lpconst', D1.lpconst, lpc, 1e-13);
assert(isequal(pdmat_inv(Phi), D1.invPhi));

% verify statistics

Mean0 = Phim / (deg - d - 1);
Mode0 = Phim / (deg + d + 1);

devcheck('mean (D0)', pdmat_fullform(mean(D0)), Mean0, 1e-13);
devcheck('mean (D1)', pdmat_fullform(mean(D1)), Mean0, 1e-13);
devcheck('mode (D0)', pdmat_fullform(mode(D0)), Mode0, 1e-13);
devcheck('mode (D1)', pdmat_fullform(mode(D1)), Mode0, 1e-13);

% verify logpdf evaluation

n = 100;
C_s = randpdm('s', d, n);
C_d = randpdm('d', d, n);
C_f = randpdm('f', d, n);

L1_s0 = D0.logpdf(C_s);
L1_d0 = D0.logpdf(C_d);
L1_f0 = D0.logpdf(C_f);

L1_s1 = D1.logpdf(C_s);
L1_d1 = D1.logpdf(C_d);
L1_f1 = D1.logpdf(C_f);

assert(isequal(size(L1_s0), [1, n]));
assert(isequal(size(L1_d0), [1, n]));
assert(isequal(size(L1_f0), [1, n]));

assert(isequal(L1_s0, L1_s1));
assert(isequal(L1_d0, L1_d1));
assert(isequal(L1_f0, L1_f1));

L0_s = my_logpdf(Phim, deg, C_s, lpc);
L0_d = my_logpdf(Phim, deg, C_d, lpc);
L0_f = my_logpdf(Phim, deg, C_f, lpc);

devcheck('loglik (s)', L1_s0, L0_s, 1e-12);
devcheck('loglik (d)', L1_d0, L0_d, 1e-12);
devcheck('loglik (f)', L1_f0, L0_f, 1e-12);

% verify sampling

N0 = 50;
N = 2e5;

C = D0.sample();
assert(is_pdmat(C) && C.d == d && C.n == 1);
assert(C.ty == 'f');
assert(isequal(C.v, C.v'));
assert(all(eig(C.v)) > 0);

Cs = D0.sample(N);
assert(is_pdmat(Cs) && Cs.d == d && Cs.n == N);
for i = 1 : N0
    Ci = pdmat_fullform(Cs, i);
    assert(isequal(Ci, Ci'));
    assert(all(eig(Ci)) > 0);
end

Cmean0 = pdmat_fullform(mean(D0));
Cmean = mean(Cs.v, 3);
devcheck('sample mean', Cmean, Cmean0, 5e-3);



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


function L = my_logpdf(Phim, m, C, lpc)

d = C.d;
n = C.n;
L = zeros(1, n);

for i = 1 : n
    Ci = pdmat_fullform(C, i);
    
    u1 = - (m + d + 1) / 2 * lndet(Ci);
    u2 = (-1/2) * trace(Ci \ Phim);
    L(i) = u1 + u2 - lpc;
end


function devcheck(name, x1, x2, thres)

d = max(abs(x1(:) - x2(:)));
if d > thres
    warning('test_wishartd:devcheck', ...
        '%s check warning with dev = %g', name, d); 
end


