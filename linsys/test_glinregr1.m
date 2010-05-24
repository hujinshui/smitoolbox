function test_glinregr1()
% A unit testing function for testing the correctness of glinregr1
% implementation
%

% Created by Dahua Lin, on April 15, 2010
%


%% main skeleton

n = 1000;

for d = [1 2 3 5]    
    for m = [0 1 3]        
        fprintf('On d = %d, m = %d, n= %d:\n', d, m, n);
        
        fprintf('    by_noiseless_est \n');
        do_test_by_noiseless_est(d, m, n);  
        
        fprintf('    by_compare_est \n');
        do_test_by_compare_est(d, m, n);
    end        
end

%% core test function

function do_test_by_noiseless_est(d, m, n)

if isequal(m, 0)
    m = 1;
    empty_weight = true;
else
    empty_weight = false;
end

c0 = randn(d, m);

x = rand(d, n * m);
y = zeros(1, n * m);
for k = 1 : m
    si = (k-1)*n+1 : k*n;
    y(si) = c0(:, k)' * x(:, si);
end

if empty_weight 
    W = [];
else
    W = l2mat(m, repnum(n * ones(1, m)));
end

glr = glinregr1(d);
ce = glr.estimate({x, y, 1}, W);

devcheck('noiseless_est_accuracy', c0, ce, 1e-8);


function do_test_by_compare_est(d, m, n)

if isequal(m, 0)
    m = 1;
    empty_weight = true;
else
    empty_weight = false;
end

c0 = rand(d, m);
x = rand(d, n * m);
y = zeros(1, n * m);
sigma = rand(1, n * m) + 0.2;
s2 = sigma.^2;

for k = 1 : m
    si = (k-1) * n+1 : k*n;
    y(si) = c0(:,k)' * x(:,si) + sigma(si) .* randn(1, n);    
end

if empty_weight
    W = [];
    W0 = ones(1, n * m);
else
    W0 = l2mat(m, repnum(n * ones(1, m)));
    W0 = W0 + rand(size(W0)) * 0.5;
    W = W0;
end

c1 = safe_estimate(x, y, W0, s2);

glr = glinregr1(d);
ce = glr.estimate({x, y, s2}, W);

devcheck('compare_est', c1, ce, 1e-12);


%% Auxiliary functions

function devcheck(title, X1, X2, thres)

dev = max(abs(X1(:) - X2(:)));
if dev > thres
    warning('test_glinregr1:largedev', ...
        'Large deviation on %s with dev = %g.', title, dev);
end


function c = safe_estimate(X, y, W, s2)

d = size(X, 1);
m = size(W, 1);

c = zeros(d, m);
for k = 1 : m
    w = W(k,:) ./ s2;
    H = bsxfun(@times, w, X) * X';
    f = X * (w .* y)';
    c(:, k) = H \ f;
end





