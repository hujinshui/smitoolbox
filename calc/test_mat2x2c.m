function test_mat2x2c(varargin)
% Test the correctness of implementation of 2 x 2 matric computation
%
%   test_mat2x2c item1 item2 ...
%       
%   item is the name of the computation to be tested, which can be
%   either of the following values:
%   'inv':      matrix inverse
%   'det':      matrix determinant
%   'chol':     Cholesky decomposition
%   'eigs':     eigen-system analysis
%   'sqrtm':    square root of matrix
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

assert(iscellstr(varargin) && ... 
    all(ismember(varargin, {'inv', 'det', 'chol', 'eigs', 'sqrtm'})), ...
    'test_mat2x2c:invalidarg', ...
    'Invalid item.');

%% main skeleton

items = varargin;
n = length(items);

for i = 1 : n
    
    item = items{i};
    
    feval(['test_' item]);
end

%% Test functions

function test_det() %#ok<DEFNU>

disp('Testing det ...');

% Test generic matrix inversion

A = gen_rand_gm_ns();
n = size(A, 3);

R1 = zeros(1, n);
for i = 1 : n
    R1(i) = det2x2(A(:,:,i));
end

R = det2x2(A);
assert(isequal(R1, R));

R0 = arrayfun(@(i) det(A(:,:,i)), 1:n);
devs = abs(R - R0);

if any(devs > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'det2x2 on gm has large deviation with max(dev) = %g.', max(devs));
end

% Test symmetric matrix inversion

B = gen_rand_sm();
nb = size(B, 2);

S1 = zeros(1, nb);
for i = 1 : nb    
    S1(i) = det2x2(B(:,i));
end

S = det2x2(B);
assert(isequal(S1, S));

B = reshape(B([1 2 2 3], :), [2 2 nb]);
S0 = arrayfun(@(i) det(B(:,:,i)), 1:nb);
devs_b = abs(S - S0);

if any(devs_b > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'det2x2 on sm has large deviation with max(dev) = %g.', max(devs_b));
end



function test_inv() %#ok<DEFNU>

disp('Testing inv ...');

% Test generic matrix inversion

A = gen_rand_gm_ns();
n = size(A, 3);

R1 = zeros(size(A));
for i = 1 : n
    R1(:,:,i) = inv2x2(A(:,:,i));
end

R = inv2x2(A);

assert(isequal(R1, R));

I = eye(2);

devs = arrayfun(@(i) calcdev(A(:,:,i) * R(:,:,i), I), 1:n) + ...
    arrayfun(@(i) calcdev(R(:,:,i) * A(:,:,i), I), 1:n);

if any(devs > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'inv2x2 on gm has large deviation with max(dev) = %g.', max(devs));
end

% Test symmetric matrix inversion

B = gen_rand_sm();
nb = size(B, 2);

S1 = zeros(2, 2, nb);
for i = 1 : nb
    r = inv2x2(B(:,i));
    S1(:,:,i) = [r(1) r(2); r(2) r(3)];
end

S = inv2x2(B);
S = reshape(S([1 2 2 3], :), [2 2 nb]);
assert(isequal(S1, S));

B = reshape(B([1 2 2 3], :), [2 2 nb]);

devs_b = arrayfun(@(i) calcdev(B(:,:,i) * S(:,:,i), I), 1:nb) + ...
    arrayfun(@(i) calcdev(S(:,:,i) * B(:,:,i), I), 1:nb);

if any(devs_b > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'inv2x2 on sm has large deviation with max(dev) = %g.', max(devs_b));
end


function test_chol()  %#ok<DEFNU>

disp('Testing chol ...');

A = gen_rand_sm();
n = size(A, 2);
A = reshape(A([1 2 2 3], :), [2 2 n]);

R1 = zeros(2, 2, n);
for i = 1 : n
    R1(:,:,i) = chol2x2(A(:,:,i));
end
R = chol2x2(A);

assert(isequal(R1, R));

B = zeros(2, 2, n);
for i = 1 : n
    B(:,:,i) = R(:,:,i) * R(:,:,i)';
end

devs = arrayfun(@(i) calcdev(A(:,:,i), B(:,:,i)), 1:n);

if any(devs > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'chol2x2 on sm has large deviation with max(dev) = %g.', max(devs));
end


function test_eigs()  %#ok<DEFNU>

disp('Testing eigs ...');

A = gen_rand_sm();
n = size(A, 2);
A = reshape(A([1 2 2 3], :), [2 2 n]);

e1 = zeros(2, n);
t1 = zeros(1, n);
for i = 1 : n
    [e1(:,i), t1(i)] = eigs2x2(A(:,:,i));
end
[e, t] = eigs2x2(A);

assert(isequal(e1, e) && isequal(t1, t));
assert(all(e(1,:) >= e(2,:)));
assert(all(abs(t) <= pi/2));

R = rotmat2(t);

B = zeros(2, 2, n);
for i = 1 : n
    B(:,:,i) = R(:,:,i)' * [e(1,i) 0; 0 e(2,i)] * R(:,:,i);
end

devs = arrayfun(@(i) calcdev(A(:,:,i), B(:,:,i)), 1:n);

if any(devs > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'eigs2x2 on sm has large deviation with max(dev) = %g.', max(devs));
end


function test_sqrtm()  %#ok<DEFNU>

disp('Testing sqrtm ...');

A = gen_rand_sm();
n = size(A, 2);
A = reshape(A([1 2 2 3], :), [2 2 n]);

R1 = zeros(2, 2, n);
for i = 1 : n
    R1(:,:,i) = sqrtm2x2(A(:,:,i));
end
R = sqrtm2x2(A);

assert(isequal(R1, R));

B = zeros(2, 2, n);
for i = 1 : n
    B(:,:,i) = R(:,:,i) * R(:,:,i);
end

devs = arrayfun(@(i) calcdev(A(:,:,i), B(:,:,i)), 1:n);

if any(devs > 1e-13)
    warning('test_mat2x2c:largedev', ...
        'sqrtm2x2 on sm has large deviation with max(dev) = %g.', max(devs));
end





%% Auxiliary functions

function A = gen_rand_gm()
% generate generic matrices

v = -1 : 0.2 : 1;
V = v(cartprod([11 11 11 11]));

A = reshape(V, 2, 2, size(V, 2));


function A = gen_rand_gm_ns()
% generate generic non-singular matrices

v = -1 : 0.2 : 1;
V = v(cartprod([11 11 11 11]));

dv = V(1,:) .* V(4,:) - V(2,:) .* V(3,:);

V = V(:, abs(dv) > 1e-8);
A = reshape(V, 2, 2, size(V, 2));


function A = gen_rand_sm()
% generate symmetric (p.d.) matrices

p = 0.1 : 0.1 : 1.5;
q = 0.1 : 0.1 : 1.5;
t = (0 : 1/16 : 1) * (2 * pi);

cp = cartprod([length(p), length(q), length(t)]);

ps = p(cp(1,:));
qs = q(cp(2,:));
ts = t(cp(3,:));

as = ps .* cos(ts).^2 + qs .* sin(ts).^2;
bs = (qs - ps) .* cos(ts) .* sin(ts);
cs = qs .* cos(ts).^2 + ps .* sin(ts).^2;

A = [as; bs; cs];
A(abs(A) < 1e-14) = 0;


function r = calcdev(X, Y)

r = max(abs(X(:) - Y(:)));

