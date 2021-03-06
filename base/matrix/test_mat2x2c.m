function test_mat2x2c(varargin)
% Test the correctness of implementation of 2 x 2 matric computation
%
%   test_mat2x2c item1 item2 ...
%       
%   item is the name of the computation to be tested, which can be
%   either of the following values:
%   'inv':      matrix inverse
%   'det':      matrix determinant
%   'trace':    matrix trace
%   'polarm':   polar representation
%   'sqrtm':    matrix square root
%   'chol':     Cholesky decomposition
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

assert(iscellstr(varargin) && ... 
    all(ismember(varargin, {'inv', 'det', 'trace', 'polarm', 'sqrtm', 'chol'})), ...
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

function test_inv() %#ok<DEFNU>

n = 1000;

X1 = gen_gm(n);

% test 4 x n form

R1 = inv2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(R1), [4 n1]));

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
R2 = inv2x2(X2);
assert(isequal(size(R2), [2 2 n1]));

assert(isequal(R1, reshape(R2, [4 n1])));

% test 3 x n form

X3 = gen_pdm(n);
n3 = size(X3, 2);
R3 = inv2x2(X3);
assert(isequal(size(R3), [3 n3]));

X4 = reshape(X3([1 2 2 3], :), [2 2 n3]);
R4 = reshape(R3([1 2 2 3], :), [2 2 n3]);

% compare result

Z20 = [ones(1, n1); zeros(1, n1); zeros(1, n1); ones(1, n1)];
Z2 = zeros(4, n1);
for i = 1 : n1
    M = 0.5 * (R2(:,:,i) * X2(:,:,i) + X2(:,:,i) * R2(:,:,i));
    Z2(:,i) = M(:);
end

Z40 = [ones(1, n3); zeros(1, n3); zeros(1, n3); ones(1, n3)];
Z4 = zeros(4, n3);
for i = 1 : n3
    M = 0.5 * (R4(:,:,i) * X4(:,:,i) + X4(:,:,i) * R4(:,:,i));
    Z4(:,i) = M(:);
end


compare('inv_gm', Z20, Z2, 5e-12);
compare('inv_sm', Z40, Z4, 5e-12);


function test_det() %#ok<DEFNU>

n = 1000;

X1 = gen_gm(n);

% test 4 x n form

r1 = det2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(r1), [1 n1]));

r10 = X1(1,:) .* X1(4,:) - X1(2,:) .* X1(3,:);

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
r2 = det2x2(X2);
assert(isequal(size(r2), [1 n1]));

assert(isequal(r1, r2));

% test 3 x n form

X3 = gen_pdm(n);
n3 = size(X3, 2);
r3 = det2x2(X3);
assert(isequal(size(r3), [1 n3]));

r30 = X3(1,:) .* X3(3,:) - X3(2,:) .* X3(2,:);


% compare result

compare('det_gm', r10, r1, 1e-15);
compare('det_sm', r30, r3, 1e-15);


function test_trace() %#ok<DEFNU>

n = 1000;

X1 = gen_gm(n);

% test 4 x n form

r1 = trace2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(r1), [1 n1]));

r10 = X1(1,:) + X1(4,:);

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
r2 = trace2x2(X2);
assert(isequal(size(r2), [1 n1]));

assert(isequal(r1, r2));

% test 3 x n form

X3 = gen_pdm(n);
n3 = size(X3, 2);
r3 = trace2x2(X3);
assert(isequal(size(r3), [1 n3]));

r30 = X3(1,:) + X3(3,:);


% compare result

compare('trace_gm', r10, r1, 1e-15);
compare('trace_sm', r30, r3, 1e-15);



function test_polarm() %#ok<DEFNU>

n = 1000;

X0 = gen_pdm(n);

% test 4 x n form

X1 = X0([1 2 2 3], :);
R1 = polarm2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(R1), [3 n1]));

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
R2 = polarm2x2(X2);
assert(isequal(size(R2), [3 n1]));

% test 3 x n form

X3 = X0;
R3 = polarm2x2(X3);
assert(isequal(size(R3), [3 n1]));

assert(isequal(R1, R2, R3));
assert(all(R1(1,:) >= R1(2,:)));


% generate ground-truth

M0 = X2;
M1 = zeros(size(M0));
for i = 1 : n1
    M1(:,:,i) = from_polar(R1(:,i));
end

% compare result

compare('polarm_sm', M0, M1, 2e-14);




function test_sqrtm() %#ok<DEFNU>

n = 1000;

X0 = gen_pdm(n);

% test 4 x n form

X1 = X0([1 2 2 3], :);
R1 = sqrtm2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(R1), [4 n1]));

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
R2 = sqrtm2x2(X2);
assert(isequal(size(R2), [2 2 n1]));

% test 3 x n form

X3 = X0;
R3 = sqrtm2x2(X3);
assert(isequal(size(R3), [3 n1]));

R2a = reshape(R2, 4, n1);
assert(isequal(R1, R2a, R3([1 2 2 3], :)));

% generate ground-truth

M0 = X2;
M1 = zeros(size(M0));
for i = 1 : n1
    M1(:,:,i) = R2(:,:,i) * R2(:,:,i);
end

% compare result

compare('sqrtm_sm', M0, M1, 2e-14);



function test_chol() %#ok<DEFNU>

n = 1000;

X0 = gen_pdm(n);

% test 4 x n form

X1 = X0([1 2 2 3], :);
R1 = chol2x2(X1);
n1 = size(X1, 2);
assert(isequal(size(R1), [4 n1]));
assert(all(R1(3, :) == 0));

% test 2 x 2 x n form

X2 = reshape(X1, [2 2 n1]);
R2 = chol2x2(X2);
assert(isequal(size(R2), [2 2 n1]));

% test 3 x n form

X3 = X0;
R3 = chol2x2(X3);
assert(isequal(size(R3), [3 n1]));

R2a = reshape(R2, 4, n1);
R3a = [R3(1:2, :); zeros(1, n1); R3(3,:)];
assert(isequal(R1, R2a, R3a));

% generate ground-truth

M0 = X2;
M1 = zeros(size(M0));
for i = 1 : n1
    M1(:,:,i) = R2(:,:,i) * R2(:,:,i)';
end

% compare result

compare('chol_sm', M0, M1, 2e-14);




%% Data generation


function X = gen_gm(n)
% generate a batch of generic matrix

X = randn(4, n);

dv = X(1,:) .* X(4,:) - X(2,:) .* X(3,:);

si = abs(dv) > 1e-3;
X = X(:, si);


function X = gen_pdm(n)
% generate a batch of positive definite matrix

X = randn(2, 2, n);

for i = 1 : n
    X(:,:,i) = X(:,:,i) * X(:,:,i)';
end

X = reshape(X, 4, n);
X = X([1 2 4], :);

dv = X(1,:) .* X(3,:) - X(2,:).^2;
si = abs(dv) > 1e-3;

X = X(:, si);


function M = from_polar(p)

a = p(1);
b = p(2);
t = p(3);

c = cos(t);
s = sin(t);

R = [c -s; s c];
M = R * [a 0; 0 b] * R';




%% Auxiliary functions

function compare(title, A, B, thres)

d = max(abs(A(:) - B(:)));
if d > thres
    warning('test_mat2x2c:largedev', ...
        'Large observation observed for %s (dev = %g)', title, d);
end


