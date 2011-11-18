function [x, lambda] = qplec(H, f, A, b)
% Solves quadratic programming problems with linear equation constraints      
%
%   x = qplec(H, f, A, b);
%       solves the following constrained quadratic programming problem:
%       
%       minimize (1/2) * x' * H * x + f' * x
%           s.t. A * x = b
%
%       Suppose there are m constraints, then size of H should be n x n,
%       and the size of A should be m x n.
%       The size of f and b can be configured in either of the following 
%       ways:
%       
%       (1) f is an n x 1 vector, and b is an m x 1 vector, specifying
%           a qp problem. In the output, the size of x is n x 1.
%
%       (2) f is an n x k matrix, and b is an m x 1 vector, specifying
%           k different qp problems, corresponding to different f.
%           In the output, the size of x is n x k, where x(:,i) 
%           corresponds to f(:, i).
%
%       (3) f is an n x 1 vector, and b is an m x k matrix, specifying
%           k different qp problems, corresponding to different b.
%           In the output, the size of x is n x k, where x(:,i)
%           corresponds to b(:, i).
%
%       (4) f is an n x k matrix, and b is an m x k matrix, specifying
%           k different qp problems, with different f and b.
%           In the output, the size of x is n x k, where x(:,i)
%           corresponds to f(:,i) and b(:,i).
%
%   [x, lambda] = qplec(H, f, A, b);
%       also returns the dual solution (lambda). 
%
%   Remarks:
%       - The function directly computes the solution using closed form 
%         solution.
%

% History
% -------
%   - Created by Dahua Lin, on Nov 26, 2009
%   - Modified by Dahua Lin, on Jul 21, 2010
%       - change the error handling to light weighting
%   - Modified by Dahua Lin, on April 16, 2011
%

%% parse and verify input arguments

if ~(isfloat(H) && isreal(H) && ndims(H) == 2 && size(H,1) == size(H,2))
    error('qplec:invalidarg', 'H should be a real valued square matrix.');
end

if ~(isfloat(f) && isreal(f) && ndims(f) == 2)
    error('qplec:invalidarg', 'f should be a real valued matrix.');
end

n = size(H, 1);
if size(f, 1) ~= n
    error('qplec:invalidarg', 'The sizes of H and f are inconsistent.');
end
kf = size(f, 2);


if ~(isfloat(A) && isreal(A) && ndims(A) == 2)
    error('qplec:invalidarg', 'A should be a real valued matrix.');
end

if ~(isfloat(b) && isreal(b) && ndims(b) == 2)
    error('qplec:invalidarg', 'b should be a real valued matrix.');
end

if size(A, 2) ~= n
    error('qplec:invalidarg', 'A should have n columns.');
end

m = size(A, 1);
if size(b, 1) ~= m
    error('qplec:invalidarg', 'b should have m columns.');
end

kb = size(b, 2);

if kf > 1 && kb > 1
    if kf ~= kb
        error('qplec:invalidarg', ...
            'when both f and b have multiple columns, the number of columns should be the same.');
    end
end



%% main

% solve the unconstraint solution

x0 = - (H \ f);

% solve the dual problem

if kf == kb
    dif = A * x0 - b;
else
    dif = bsxfun(@minus, A * x0, b);
end

G = H \ A';
Q = A * G;
Q = (Q + Q') / 2;
lambda = Q \ dif;

% derive the primal solution

if size(x0, 2) == size(lambda, 2)
    x = x0 - G * lambda;
else
    x = bsxfun(@minus, x0, G * lambda);
end

