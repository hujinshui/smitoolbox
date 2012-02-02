function f = glinreg_objfun(X, Y, w, rho, rc)
%GLINREG_OBJFUN Generalized linear regression objective function
%
%   A generalized linear regression problem is an optimization problem
%   formulated as minimizing the following objective
%
%       f(theta, theta0) = 
%           \sum_{i=1}^n + w_i * rho((theta' * x_i + theta_0) - y_i)
%           + (rc/2) * ||theta||^2
%
%   Suppose the dimension of the feature space is d, and that of the
%   output space is q. Then each x_i is a d x 1 vector, and each y_i
%   is a q x 1 vector, and theta is a d x q matrix, and theta_0 is 
%   a q x 1 vector.
%
%   Here, rho is a function that takes a q x 1 difference vector as
%   input and yields a loss value.
%   
%   The solution is a vector of length (d+1) x q. Suppose sol is 
%   the solution vector, then
%
%   if q == 1
%       theta = sol(1:d);
%       theta0 = sol(d+1);
%   else
%       A = reshape(sol, d+1, q);
%       theta = A(1:d, :);
%       theta0 = A(d+1, :);
%   end
%
%
%   f = GLINREG_OBJFUN(X, Y, w, rho, rc);
%   
%       Returns the objective function f for the generalized linear
%       regression problem formulated above.
%
%       Input arguments:
%
%       - X:        The feature matrix, size: d x n
%       
%       - Y:        The matrix of expected output, size: q x n
%
%       - w:        The sample weights, empty or a 1 x n vector.
%
%       - rho:      The loss function handle, which takes difference
%                   vectors as input, and outputs the loss values.
%
%       - rc:       The regularization coefficient: a scalar or a d x 1
%                   vector.
%
%       Output arguments:
%       - f:        The objective function handle, which can be used as
%
%                   v = f(sol);
%                   [v, g] = f(sol);
%
%                   Here, v and g are respectively the functon value and
%                   gradient evaluated at sol.
%

% Created by Dahua Lin, on Jan 15, 2012
%

%% verify input

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('glinreg_objfun:invalidarg', 'X should be a real matrix.');
end

if ~(isfloat(Y) && isreal(Y) && ndims(Y) == 2)
    error('glinreg_objfun:invalidarg', 'X should be a real matrix.');
end

n = size(X, 2);
if size(Y, 2) ~= n
    error('glinreg_objfun:invalidarg', 'X and Y have different number of columns.');
end

if ~isempty(w)
    if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('glinreg_objfun:invalidarg', 'w should be a real vector of length n.');
    end
end

if ~isa(rho, 'function_handle')
    error('glinreg_objfun:invalidarg', 'rho should be a function handle.');
end

if ~(isfloat(rc) && isscalar(rc) && isreal(rc))
    error('glinreg_objfun:invalidarg', 'rc should be a real scalar.');
end


%% main

n = size(X, 2);
Xa = [X; ones(1, n)];
f_loss = comb_lossfun(Xa, Y, size(Y, 1),  w, @(z, y) rho(z - y));

if rc == 0
    f = f_loss;
else
    d = size(X, 1);
    q = size(Y, 1);
    c = [ones(d, 1); 0];
    if q > 1
        c = c(:, ones(1,q));
        c = c(:);
    end 
    f = comb_objfun(1, f_loss, rc, tikregf(c));
end

