function test_pdmat()
% Test the correctness of the implementation of pdmat_* functions
%
%   test_pdmat;
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% main

tys = {'s', 'd', 'f'};
ds = [1, 2, 5];
ns = [1, 3];

test_construction(tys, ds, ns);

run_tests('Test pick & fullform', tys, ds, ns, @test_pick_and_fullform);
run_tests('Test diag',  tys, ds, ns, @test_diag);
run_tests('Test scale', tys, ds, ns, @test_scale);
run_tests('Test inv', tys, ds, ns, @test_inv);
run_tests('Test lndet', tys, ds, ns, @test_lndet);

run_tests('Test mvmul', tys, ds, ns, @test_mvmul);
run_tests('Test lsolve', tys, ds, ns, @test_lsolve);
run_tests('Test quad', tys, ds, ns, @test_quad);
run_tests('Test pwquad', tys, ds, 1, @test_pwquad);
run_tests('Test choltrans', tys, ds, 1, @test_choltrans);

test_plus(tys, ds, ns);

% test_plus(tys, ds, ns);

disp('All test cases pass!');
disp(' ');


%% test cases

function test_construction(tys, ds, ns)

disp('Test construction ...');

% test construction with full spec
for t = 1 : numel(tys)
    ty = tys{t};
    for d = ds
        for n = ns
            randpdm(ty, d, n);
        end
    end
end

% test construction with simplified spec

pdm_s = randpdm('s', 1, 1);
assert(isequal(pdmat(pdm_s.v), pdm_s));

pdm_d = randpdm('d', 2, 1);
assert(isequal(pdmat(pdm_d.v), pdm_d));

pdm_f = randpdm('f', 2, 1);
assert(isequal(pdmat(pdm_f.v), pdm_f));


function test_pick_and_fullform(ty, d, n)

S = randpdm(ty, d, n);
for i = 1 : n
    Si = pdmat_pick(S, i);
    check_valid(Si, ty, d, 1);
    
    Fi0 = pdmat_fullform(S, i);
    Fi1 = pdmat_fullform(Si);
    assert(isa(Fi0, 'double') && isequal(size(Fi0), [d d]));
    assert(isa(Fi1, 'double') && isequal(size(Fi1), [d d]));
    assert(isequal(Fi0, Fi1));
end


function test_diag(ty, d, n)

S = randpdm(ty, d, n);
V0 = zeros(d, n);
for i = 1 : n
    Ci = pdmat_fullform(S, i);
    V0(:, i) = diag(Ci);
end
V1 = pdmat_diag(S);
assert(isequal(V0, V1));


function test_scale(ty, d, n)

S = randpdm(ty, d, n);
c = pi;
R = pdmat_scale(S, c);
check_valid(R, ty, d, n);

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    Rm = pdmat_fullform(R, i);
    assert(norm(Sm * c - Rm) < 1e-15);    
end

if n == 1
    K = 5;
    c = rand(1, K);
    R = pdmat_scale(S, c);
    check_valid(R, ty, d, K);
    
    Sm = pdmat_fullform(S);
    for i = 1 : K
        Rm = pdmat_fullform(R, i);
        assert(norm(Sm * c(i) - Rm) < 1e-15);
    end
    
else
    c = rand(1, n);
    R = pdmat_scale(S, c);
    check_valid(R, ty, d, n);
    
    for i = 1 : n
        Sm = pdmat_fullform(S, i);
        Rm = pdmat_fullform(R, i);
        assert(norm(Sm * c(i) - Rm) < 1e-15);
    end    
end



function test_inv(ty, d, n)

S = randpdm(ty, d, n);
R = pdmat_inv(S);
check_valid(R, ty, d, n);

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    Rm = pdmat_fullform(R, i);
    
    assert(norm(Sm * Rm - eye(d)) < 1e-12);
end


function test_lndet(ty, d, n)

S = randpdm(ty, d, n);
r = pdmat_lndet(S);
assert(isa(r, 'double') && isreal(r));
assert(isequal(size(r), [1, n]));

for i = 1: n
    Sm = pdmat_fullform(S, i);
    cr0 = lndet(Sm);
    assert(abs(r(i) - cr0) < 1e-12);
end


function test_mvmul(ty, d, n)

S = randpdm(ty, d, n);

X = rand(d, n);
Y = pdmat_mvmul(S, X);
assert(isa(Y, 'double') && isequal(size(Y), [d, n]));

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    y0 = Sm * X(:, i);
    assert(norm(y0 - Y(:,i)) < 1e-13);    
end

if n > 1
    return;
end

nx = 5;
X = rand(d, nx);
Y = pdmat_mvmul(S, X);
assert(isa(Y, 'double') && isequal(size(Y), [d, nx]));

Y0 = pdmat_fullform(S) * X;
assert(norm(Y - Y0) < 1e-12);
    

function test_lsolve(ty, d, n)

S = randpdm(ty, d, n);

X = rand(d, n);
Y = pdmat_lsolve(S, X);
assert(isa(Y, 'double') && isequal(size(Y), [d, n]));

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    y0 = Sm \ X(:, i);
    assert(norm(y0 - Y(:,i)) < 1e-12);
end

if n > 1
    return;
end

nx = 5;
X = rand(d, nx);
Y = pdmat_lsolve(S, X);
assert(isa(Y, 'double') && isequal(size(Y), [d, nx]));

Y0 = pdmat_fullform(S) \ X;
assert(norm(Y - Y0) < 1e-12);


function test_quad(ty, d, n)

S = randpdm(ty, d, n);

x = rand(d, 1);
y = rand(d, 1);

rx = pdmat_quad(S, x);
rxy = pdmat_quad(S, x, y);

assert(isa(rx, 'double') && isequal(size(rx), [n, 1]));
assert(isa(rxy, 'double') && isequal(size(rxy), [n, 1]));

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    
    crx0 = x' * Sm * x;
    crxy0 = x' * Sm * y;
    
    assert(abs(rx(i) - crx0) < 1e-13);
    assert(abs(rxy(i) - crxy0) < 1e-13);
end

nx = 5;
X = rand(d, nx);
Y = rand(d, nx);

Rxx = pdmat_quad(S, X);
Rxy = pdmat_quad(S, X, Y);

assert(isa(Rxx, 'double') && isequal(size(Rxx), [n, nx]));
assert(isa(Rxy, 'double') && isequal(size(Rxy), [n, nx]));

for i = 1 : n
    Sm = pdmat_fullform(S, i);
    
    crx0 = dot(X, Sm * X, 1);
    crxy0 = dot(X, Sm * Y, 1);
    
    assert(norm(Rxx(i, :) - crx0) < 1e-12);
    assert(norm(Rxy(i, :) - crxy0) < 1e-12);
end


function test_pwquad(ty, d, n)

assert(n == 1);
S = randpdm(ty, d, n);

nx = 5; X = rand(d, nx);
ny = 6; Y = rand(d, ny);

Rxx = pdmat_pwquad(S, X);
Rxy = pdmat_pwquad(S, X, Y);

assert(isa(Rxx, 'double') && isequal(size(Rxx), [nx, nx]));
assert(isa(Rxy, 'double') && isequal(size(Rxy), [nx, ny]));

Sm = pdmat_fullform(S);
Rxx0 = X' * Sm * X;
Rxy0 = X' * Sm * Y;

assert(norm(Rxx - Rxx0) < 1e-12);
assert(norm(Rxy - Rxy0) < 1e-12);


function test_choltrans(ty, d, n)

assert(n == 1);
S = randpdm(ty, d, n);

nx = 5; X = rand(d, nx);

R = pdmat_choltrans(S, X);
assert(isa(R, 'double') && isequal(size(R), [d, nx]));

Sm = pdmat_fullform(S);
R0 = chol(Sm, 'lower') * X;

assert(norm(R - R0) < 1e-12);


function test_plus(tys, ds, ns)

disp('Test plus ...');

for t1 = 1 : numel(tys)
    for t2 = 1 : numel(tys)
        
        ty1 = tys{t1};
        ty2 = tys{t2};
        
        for d = ds
            for n = ns
                
                if n == 1
                    test_plus_a(d, ty1, 1, ty2, 1);
                else
                    test_plus_a(d, ty1, n, ty2, n);
                    test_plus_a(d, ty1, 1, ty2, n);
                    test_plus_a(d, ty1, n, ty2, 1);
                end
                
            end
        end
    end    
end


function test_plus_a(d, ty1, n1, ty2, n2)

S1 = randpdm(ty1, d, n1);
S2 = randpdm(ty2, d, n2);

c1 = 2;
c2 = 3;

R = pdmat_plus(S1, S2);
C = pdmat_plus(S1, S2, c1, c2);

nr = max(n1, n2);
tyr = max_form_ty(ty1, ty2);

check_valid(R, tyr, d, nr);
check_valid(C, tyr, d, nr);

for i = 1 : nr
    if n1 == 1
        Sm1 = pdmat_fullform(S1);
    else
        Sm1 = pdmat_fullform(S1, i);
    end
    
    if n2 == 1
        Sm2 = pdmat_fullform(S2);
    else
        Sm2 = pdmat_fullform(S2, i);
    end
    
    Rm = pdmat_fullform(R, i);
    Rm0 = Sm1 + Sm2;
    assert(norm(Rm - Rm0) < 1e-13);
    
    Cm = pdmat_fullform(C, i);
    Cm0 = c1 * Sm1 + c2 * Sm2;
    assert(norm(Cm - Cm0) < 1e-13);
end



%% auxiliary functions

function S = randpdm(ty, d, n)
% generate a random pdmat struct and verify its validity

mg = 0.2; % margin to keep away from zeros

switch ty
    case 's'
        v = rand(1, n) + mg;
    case 'd'
        v = rand(d, n) + mg;
    case 'f'        
        v = zeros(d, d, n);
        for i = 1 : n
            r = orth(rand(d, d));
            cv = r' * bsxfun(@times, rand(d,1) + mg, r);
            cv = 0.5 * (cv + cv');
            v(:,:,i) = cv;
        end
    otherwise
        error('test_pdmat:rterror', 'Unknown form type %s', ty);
end

S = pdmat(ty, d, v);

% verify
check_valid(S, ty, d, n);
assert(isa(S.v, 'double') && isequal(S.v, v));


function check_valid(S, ty, d, n)

assert(is_pdmat(S));
assert(ischar(S.tag) && strcmp(S.tag, 'pdmat'));
assert(ischar(S.ty) && strcmp(S.ty, ty));
assert(isa(S.d, 'double') && isequal(S.d, d));
assert(isa(S.n, 'double') && isequal(S.n, n));

switch ty
    case 's'
        assert(isequal(size(S.v), [1, S.n]));
    case 'd'
        assert(isequal(size(S.v), [S.d, S.n]));
    case 'f'
        if S.n == 1
            assert(isequal(size(S.v), [S.d, S.d]));
        else
            assert(isequal(size(S.v), [S.d, S.d, S.n]));
        end
end


function run_tests(title, tys, ds, ns, tfunc)

fprintf('%s ...\n', title);

for t = 1 : numel(tys)
    ty = tys{t};
    for d = ds
        for n = ns
            tfunc(ty, d, n);
        end
    end
end


function tyr = max_form_ty(ty1, ty2)

s = [ty1, ty2];

switch s
    case 'ss'
        tyr = 's';
    case {'sd', 'ds', 'dd'}
        tyr = 'd';
    case {'sf', 'df', 'fs', 'fd', 'ff'}
        tyr = 'f';
    otherwise
        error('test_pdmat:rterror', 'Invalid types.');
end



