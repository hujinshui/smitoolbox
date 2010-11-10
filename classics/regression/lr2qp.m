function [H, f] = lr2qp(A, y, imat, r)
% Constructs the QP program for linear regression
%
%   [H, f] = lrqp(A, y);
%       converts the following least square linear regression problem
%
%       minimize (1/2) * ||Ax - y||^2
%
%       to the quadratic programming problem as
%
%       minimize (1/2) * x' * H * x + f' * x.
%
%       H and f are output. Here, A and y should have the same number of rows.
%       Suppose A is an m x n matrix, and y is an m x k matrix, then 
%       H is an n x n matrix, and f is an n x k matrix. In particular,
%       f(:,i) corresponds to y(:,i).
%
%   [H, f] = lrqp(A, y, imat);
%       converts the following problem (with modified norm) to QP.
%
%       minimize (1/2) * (Ax - y)' * imat * (Ax - y).
%
%       Here imat is an m x m semi-definite matrix. It can be specified
%       as either an m x m matrix, an m x 1 vector (for diagonal matrix)
%       or a scalar.
%
%       Omitting imat is equivalent to setting it to be 1.
%
%   [H, f] = lrqp(A, y, imat, r);
%       converts the following ridge regression problem to QP.
%
%       minimize (1/2) (Ax - y)' * imat * (Ax - y) + .(1/2) * x' * r * x
%
%       Here, r is the Tikhonov regularization coefficient matrix.
%       It is an n x n matrix, which can also be specified as an n x 1
%       vector (for diagonal matrix), or a scalar.
%

%  Created by Dahua Lin, on Nov 26, 2009
%

%% parse and verify input arguments

assert(isfloat(A) && ndims(A) == 2, 'linrd:invalidarg', ...
    'A should be a 2D numeric matrix.');

[m, n] = size(A);

assert(isfloat(y) && ndims(y) == 2 && size(y,1) == m, ...
    'linrd:invalidarg', 'y should be a numeric matrix with same number of rows as A.');

if nargin < 3
    mf = 0;
    imat = 1;    
    
else
    assert(isfloat(imat) && isreal(imat) && ndims(imat) == 2, ...
        'linrd:invalidarg', 'imat should be a real numeric matrix');
    
    if isscalar(imat)
        mf = 0;
        
    elseif size(imat, 1) == m && size(imat, 2) == 1
        mf = 1;
        
    elseif size(imat, 1) == m && size(imat, 2) == m
        mf = 2;
        
    else
        error('linrd:invalidarg', 'The size of imat is invalid.');
    end
end


if nargin < 4
    r = 0;
    rf = 0;
    
else
    assert(isfloat(r) && isreal(r) && ndims(r) == 2, ...
        'linrd:invalidarg', 'r should be a real numeric matrix.');
    
    if isscalar(r) 
        rf = 0;
        
    elseif size(r, 1) == n && size(r, 2) == 1
        rf = 1;
        
    elseif size(r, 1) == n && size(r, 2) == n
        rf = 2;
        
    else
        error('linrd:invalidarg', 'The size of r is invalid.');
    end
end


        
%% main

% process data term

if mf == 0
    if imat == 1
        B = A';
    else
        B = A' * imat;
    end
    
elseif mf == 1
    B = bsxfun(@times, A, imat)';
    
else
    B = A' * imat;
    
end

H = B * A;
f = -(B * y);


% regularize

if rf == 0
    if r > 0
        H(1:(n+1):n*n) = H(1:(n+1):n*n) + r;
    end
elseif rf == 1    
    H(1:(n+1):n*n) = H(1:(n+1):n*n) + r';
else
    H = H + r;
end
    
   
