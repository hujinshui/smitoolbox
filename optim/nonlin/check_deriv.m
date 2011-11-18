function [r1, r2] = check_deriv(f, X, h)
% A function for checking whether derivatives are computed correctly
%
%   r1 = check_deriv(f, X);
%   r1 = check_deriv(f, X, h);
%       
%       This statement takes an objective function f and a set of vectors
%       in X as input. Here, f should be a function handle that supports
%       the following usage:
%
%           [v, g] = f(x);
%
%       Here, f returns the objective value and the gradient at x.
%
%       Supppose the input dimension for f is d, and there are n check 
%       points. Then X should be a matrix of size d x n, with each column
%       being a test vector.
%
%       This function will compute the gradient of f at these test vectors
%       and compare them with those yielded by finite difference
%       approximation. In addition, one can specify h, which is the
%       difference value used in approximation. By default, it is 1e-5.
%
%       In the output, r1 is a vector of size 1 x n. In particular, r1(i)
%       is the norm of the difference between the gradient yielded by f
%       and that by approximation.
%
%   [r1, r2] = check_deriv(f, X);
%   [r1, r2] = check_deriv(f, X, h);
%
%       This statement also checks the Hessian matrices. Here, f should
%       support the following usage:
%
%           [v, g, H] = f(x);
%
%       Here, f returns the objective value, the gradient, as well as the
%       Hessian matrix at x.
%
%       In the output, r1 gives the norm of difference between gradients,
%       while r2 gives the norm of difference between Hessian matrices.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% verify input arguments

if ~isa(f, 'function_handle')
    error('check_deriv:invalidarg', 'f should be a function handle.');
end

if ~(isfloat(X) && ndims(X) == 2)
    error('check_deriv:invalidarg', 'X should be a numeric matrix.');
end

if nargin < 3
    h = 1e-5;
else
    if ~(isfloat(h) && isscalar(h) && h > 0)
        error('check_deriv:invalidarg', 'h should be a positive scalar.');
    end
end

%% main

n = size(X, 2);

r1 = zeros(1, n);

if nargout < 2
    for i = 1 : n
        x = X(:,i);
        [v, g] = f(x); %#ok<ASGLU>
        ga = approx_deriv(f, x, h);
        r1(i) = norm(g - ga);
    end
else
    r2 = zeros(1, n);
    for i = 1 : n
        x = X(:,i);
        [v, g, H] = f(x); %#ok<ASGLU>
        [ga, Ha] = approx_deriv(f, x, h);
        r1(i) = norm(g - ga);
        r2(i) = norm(H - Ha);
    end
end

