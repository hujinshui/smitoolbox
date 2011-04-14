function P = linear_svm_prob(X, y, c, op)
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
%   This function also supports a LP-variation for large-scale application,
%   where the objective function to be minimize is modified to be
%
%       ||w||_1 + \sum_i C xi_i
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
%   P = linear_svm_prob(X, y, c, 'Lp');
%
%       constructs the variant LP-problem for linear SVM that uses L1-norm
%       of w.
%
%       The output P is a lp_problem struct.
%
%       The solution to the constructed problem is of size 2d + 1 + n.
%       Let sol be the solution vector:
%       sol(1:d)     -  w
%       sol(d+1:2*d) -  abs(w)
%       sol(2*d+1)   -  b
%       sol(2*d+2:2*d+n+1) - xi (slack variables)
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

if nargin >= 4
    if ~(ischar(op) && strcmpi(op, 'Lp'))
        error('linear_svm_prob:invalidarg', ...
            'The 4th argument must be ''Lp'' if given.');
    end
    use_Lp = true;
else
    use_Lp = false;
end


%% main

if ~use_Lp
    
    H = sparse(1:d, 1:d, 1, d+1+n, d+1+n);
    f = [zeros(d+1, 1); c * ones(n, 1)];
    
    Xy = bsxfun(@times, X, y);
    if issparse(Xy) || n > 3 * d  % use sparse construction
        A = [Xy.', y.', speye(n, n); sparse(n, d+1), speye(n, n)];
    else  % use dense construction
        A = [Xy.', y.', eye(n, n); zeros(n, d+1), eye(n, n)];
    end
    b = [ones(n, 1); zeros(n, 1)];
    
    P = qp_problem(H, f, -A, -b, [], [], [], []);

else
        
    f = [zeros(d, 1); ones(d, 1); 0; c * ones(n, 1)];
    
    Xy = bsxfun(@times, X, y);
    
    if issparse(Xy) || n > 3 * d  % use sparse construction
        A = [ ...
            Xy.', sparse(n, d), y.', speye(n, n); ...
            sparse(n, 2*d+1), speye(n, n); ...
            speye(d, d), speye(d, d), sparse(d, 1+n); ...
            - speye(d, d), speye(d, d), sparse(d, 1+n) ];        
    else
        A = [ ...
            Xy.', zeros(n, d), y.', eye(n, n); ...
            zeros(n, 2*d+1), eye(n, n); ...
            eye(d, d), eye(d, d), zeros(d, 1+n); ...
            -eye(d, d), eye(d, d), zeros(d, 1+n)];
    end
    b = [ones(n, 1); zeros(n + 2 *d, 1)];
    
    P = lp_problem(f, -A, -b, [], [], [], []);
end


