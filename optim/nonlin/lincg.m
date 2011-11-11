function [x, flag] = lincg(A, b, x0, tol, max_iter)
%LINCG Solves a linear equation using conjugate gradient
%
%   x = LINCG(A, b, x0);
%   x = LINCG(A, b, x0, tol);
%   x = LINCG(A, b, x0, tol, max_iter);
%
%       solves a linear equation A * x = b using conjugate gradient
%       method.
%
%       Input arguments:
%       - A:            the linear coefficient matrix, or it can be a
%                       function handle such that A(x) yields A * x.
%       - b:            the vector on the right hand side of the equation
%       - x0:           initial guess of the solution
%       - tol:          the tolerance of residue norm at convergence 
%                       (when omitted, tol = 1e-6 by default)
%       - max_iter:     the maximum number of iterations
%                       (when omitted, max_iter = 100 by default)
%
%       Output arguments:
%       - x:            the solution
%
%   [x, flag] = LINCG( ... );
%       additionally returns the flag about the result of the solution.
%       Here is the list of flag meanings:
%       0:  the algorithm converges to the desired tolerance
%       1:  the algorithm iterates max_iter times but have not converged.
%       2:  early termination due to numerical issues
%

% Created by Dahua Lin, on Oct 31, 2011
%


%% verify input arguments

if isnumeric(A)
    if ~(isfloat(A) && ndims(A) == 2 && isreal(A) && size(A,1) == size(A,2))
        error('lincg:invalidarg', 'A should be a real matrix.');    
    end    
    n = size(A, 1);
    is_fun = 0;
elseif isa(A, 'function_handle')
    n = 0;      % indicate it is unknown from A
    is_fun = 1;
end

if ~(isfloat(b) && ndims(b) == 2 && size(b,2) == 1 && isreal(b) && ~issparse(b))
    error('lincg:invalidarg', 'b should be a non-sparse real vector.');
end
    
if n > 0
    if size(b,1) ~= n
        error('lincg:invalidarg', 'The dimensions of A and b are inconsistent.');
    end
else
    n = size(b, 1);
end

if ~(isfloat(x0) && isequal(size(x0), [n 1]) && isreal(x0) && ~issparse(x0))
    error('lincg:invalidarg', 'x0 should be a non-sparse real vector.');
end

if nargin < 4
    tol = 1e-6;
else
    if ~(isfloat(tol) && isscalar(tol) && isreal(tol) && tol > 0)
        error('lincg:invalidarg', 'tol should be a positive real value.');
    end
end

if nargin < 5
    max_iter = 100;
else
    if ~(isnumeric(max_iter) && isscalar(max_iter) && ...
            max_iter == fix(max_iter) && max_iter >= 1)
        error('lincg:invalidarg', 'max_iter should be a positive integer.');
    end
end

%% main

x = x0;
if is_fun
    r = b - A(x);
else
    r = b - A * x;
end
p = r;
k = 0;
converged = false;

s = r' * r;

while k < max_iter
    k = k + 1;
    
    if is_fun
        z = A(p);
    else
        z = A * p;    
    end
    alpha = s / (p' * z);    
    x = x + alpha * p;
    
    r = r - alpha * z;    
    s_pre = s;
    s = r' * r;
            
    if sqrt(s) < tol
        converged = true;
        break;
    end
    
    beta = s / s_pre;
    p = r + beta * p;
end

if converged 
    flag = 0;
else
    flag = 1;
end




