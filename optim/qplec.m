function x = qplec(H, f, A, b)
% Solves quadratic programming problems with linear equation constraints
%
%   x = qplec(H, f);
%       solves the unconstrained quadratic programming problem as
%
%       minimize (1/2) * x' * H * x + f' * x
%
%       Let x be in an n-dimensional space, then H and f should respectively
%       be an n x n positive definite matrix and an n x 1 vector.
%
%       As an extension, one can also input f as an n x k matrix, specifying
%       k different qp programs which share the same H matrix. In this case,
%       x will be an n x k matrix, with x(:,i) being the solution to the 
%       problem whose linear term coefficients are given by f(:, i).       
%
%   x = qplec(H, f, A, b);
%       solves the following constrained quadratic programming problem:
%       
%       minimize (1/2) * x' * H * x + f' * x
%           s.t. A * x = b
%
%       Suppose there are m constraints, then size of H should be n x n,
%       and the size of A should be m x n.
%       The size of f and b can be configured in either of the following ways:
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
%   Remarks:
%       - The function directly computes the solution using closed form 
%         formula.
%

% History
% -------
%   - Created by Dahua Lin, on Nov 26, 2009
%

%% parse and verify input arguments

assert(isfloat(H) && isreal(H) && ndims(H) == 2 && size(H,1) == size(H,2), ...
    'qplec:invalidarg', 'H should be a real valued square matrix.');

assert(isfloat(f) && isreal(f) && ndims(f) == 2, ...
    'qplec:invalidarg', 'f should be a real valued matrix.');

n = size(H, 1);
assert(size(f, 1) == n, 'qplec:invalidarg', 'f should have n rows.');
kf = size(f, 2);

if nargin < 3 || isempty(A)
    A = [];
    b = [];
else
    assert(isfloat(A) && isreal(A) && ndims(A) == 2, ...
        'qplec:invalidarg', 'A should be a real valued matrix.');
    
    assert(isfloat(b) && isreal(b) && ndims(b) == 2, ...
        'qplec:invalidarg', 'b should be a real valued matrix.');
    
    assert(size(A, 2) == n, 'qplec:invalidarg', 'A should have n columns.');
    
    m = size(A, 1);
    assert(size(b, 1) == m, 'qplec:invalidarg', 'b should have m columns.');
    
    kb = size(b, 2);
    
    if kf > 1 && kb > 1
        assert(kf == kb, 'qplec:invalidarg', ...
            'when both f and b have multiple columns, the number of columns should be the same.');
    end
end


%% main

x = - (H \ f);

if ~isempty(A)   % correction
    
    if kf == kb
        D = A * x - b;
    else
        D = bsxfun(@minus, A * x, b);
    end
    
    if nnz(D) > 0   % solution do not satisfy constraint, need correction
        
        G = A * (H \ A');
        G = 0.5 * (G + G');
        
        lambda = G \ D;
        
        if size(f, 2) == size(lambda, 2)
            cf = f + A' * lambda;
        else
            cf = bsxfun(@plus, f, A' * lambda);
        end
        
        x = - (H \ cf);
    end        
end


