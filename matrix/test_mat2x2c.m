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







%% Auxiliary functions

function compare(title, A, B, thres)

d = max(abs(A(:) - B(:)));
if d > thres
    warning('test_mat2x2c:largedev', ...
        'Large observation observed for %s (dev = %g)', title, d);
end


