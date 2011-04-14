function P = linear_svm_prob(X, y, c)
% Constructs the primal problem for standard linear SVM
%
%   The standard linear svm problem is formulated as follows
%
%       minimize (1/2) * ||w||^2 + \sum_i C xi_i
%           s.t. y_i (w' * x_i + b) >= 1 - xi_i
%                xi_i >= 0
%
%       Here, w is the coefficient vector, and xi_i is the slack variable.
%
%
%   P = linear_svm_prob(X, y, c);
%
%       constructs the linear SVM QP problem as formulated above.
%
%       Input arguments:
%       - X:        the sample matrix of size d x n, with each column 
%                   giving a sample.
%
%       - y:        the label vector of size 1 x n, with y(i) corresponding
%                   to X(:,i).
%
%       - c:        the weight to penalize margin deviation.
%                   
%       The output P is a qp_problem struct.
%
%       The solution to the constructed problem is of size d + 1 + n.
%       Let sol be the solution vector, then sol(1:d) is the coefficient
%       vector, sol(d+1) is b, and sol(d+2:d+n+1) is the slack variables.
%

% Created by Dahua Lin, on April 14, 2011
%

%% verify input arguments

if ~(isnumeric(X) && ndims(X) == 2 && isreal(X) && ~isempty(X))
    error('linear_svm_prob:invalidarg', ...
        'X should be an non-empty real matrix.');
end
if ~isa(X, 'double'); X = double(X); end
[d, n] = size(X);

if ~(isnumeric(y) && isreal(y) && isequal(size(y), [1 n]))
    error('linear_svm_prob:invalidarg', ...
        'y should be a numeric row vector of size 1 x n.');
end
if ~isa(y, 'double'); y = double(y); end

if ~(isnumeric(c) && isscalar(c) && isreal(c) && c > 0)
    error('linear_svm_prob:invalidarg', ...
        'c should be a numeric positive real scalar.');
end
c = double(c);

%% main

ds = d + 1 + n;     % dimension of solution space
H = sparse(1:d, 1:d, 1, ds, ds);
f = [zeros(d+1, 1); c * ones(n, 1)];

Xy = bsxfun(@times, X, y);
if issparse(Xy) || n > 3 * d  % use sparse construction
    A = [Xy.', y.', speye(n, n); sparse(n, d+1), speye(n, n)];    
else  % use dense construction
    A = [Xy.', y.', eye(n, n); zeros(n, d+1), eye(n, n)];
end
b = [ones(n, 1); zeros(n, 1)];

P = qp_problem(H, f, -A, -b, [], [], [], []);


