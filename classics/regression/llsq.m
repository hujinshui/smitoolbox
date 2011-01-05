function [o1, o2] = llsq(X, y, w, r)
% Linear leasr square problem
%
%   The function solves the (weighted) linear least square problem, as
%
%    minimize (1/2) * sum_i w_i ||x_i' * a - y_i||^2 + (1/2) * r * ||a||^2
%
%   It is well known that the solution to this problem is
%
%       optimal a = inv(X' * W * X + rI) * (X' * W * y).
%
%   
%   a = llsq(X, y);
%   a = llsq(X, y, w);
%   a = llsq(X, y, w, r);
%
%       solves the problem as formalized above and returns the optimal a.
%
%       Input arguments:
%       - X:    the design matrix. Suppose there are n independent
%               variables, each with d components, then X should be
%               an n x d matrix.
%
%       - y:    the response vector/matrix. In standard formulation,
%               y should be an n x q vector, with y(i,:) corresponding to
%               X(i, :). Here, q is the vector space dimension for the
%               response.
%
%       - w:    the weights. w can be in either of the following form:
%               - a scalar: All rows in X and y have the same weight w.
%               - a vector of length n: The i-th row has weight w(i).
%               - an positive semi-definite n x n matrix: This is a 
%                 generalization of the standard formulation which allows 
%               correlation between rows.
%
%               Omitting w is equivalent to setting w to 1, in which case,
%               the function essentially solves the ordinary least square
%               problem.
%                 
%       - r:    the regularization coefficient. r can be either of the
%               following: 
%               - a scalar: as in the formulation above
%               - a vector of length d: assigning different regularization
%                 coefficients for different components
%               - a postive definite matrix of size d x d: performing 
%                 generic Tikhonov regularization, in which the
%                 regularization term is (1/2) * a' * r * a.
%
%               Omitting r is equivalent to setting r to 0, in which case,
%               no regularization is performed.
%
%   [H, f] = llsq(X, y, w);
%       
%       returns the Hessian matrix and gradient vector(s) instead of 
%       solving the problem.
%
%       Here, H = X' * W * X + rI, and f = X' * W * y.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 5, 2011
%

%% verify input

if ~(isfloat(X) && ndims(X) == 2)
    error('llsq:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 1);

if ~(isfloat(y) && ndims(y) == 2)
    error('llsq:invalidarg', 'y should be a numeric matrix.');
end

if n ~= size(y, 1)
    error('llsq:invalidarg', 'X and y should have the same number of rows.');
end

if nargin < 3
    w = 1;
else
    if ~(isfloat(w) && ndims(w) == 2)
        error('llsq:invalidarg', 'w should be a numeric scalar/vector/matrix.');
    end
end

if nargin < 4
    r = 0;
else
    if ~(isfloat(r) && ndims(r) == 2)
        error('llsq:invalidarg', 'r should be a numeric scalar/vector/matrix.');
    end
end

%% main

% compute H and f

if isscalar(w)
    H = X' * X;
    f = X' * y;
    if w ~= 1
        H = w * H;
        f = w * f;
    end
    
elseif isvector(w)
    if numel(w) ~= n
        error('llsq:invalidarg', 'The length of w is invalid.');
    end
    if size(w, 2) > 1; w = w.'; end % turn w into a column
    H = X' * bsxfun(@times, w, X);
    H = 0.5 * (H + H');
    f = X' * bsxfun(@times, w, y);
    
else 
    if ~(size(w,1) == n && size(w, 2) == n)
        error('llsq:invalidarg', 'The size of w is invalid.');
    end
    H = X' * (W * X);
    H = 0.5 * (H + H');
    f = X' * (W * y);
    
end

% regularize H

if ~isequal(r, 0)    
    if isscalar(r) || isvector(r)
        H = adddiag(H, r);
    else
        H = H + r;
    end
end
        

% solve and return

if nargout < 2
    o1 = H \ f;
else
    o1 = H;
    o2 = f;
end

