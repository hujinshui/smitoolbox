function [H, f] = lr2qp(A, y, Q, R)
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
%       In mathematics, the relations between (H, f) and (A, y) can be
%       given by the following formulas:
%
%           H = A' * A   and    f = - A' * y
% 
%       Here, A and y should be matrices with the same number of rows.
%       Suppose A is an m x n matrix, and y is an m x k matrix, then 
%       H is an n x n matrix, and f is an n x k matrix. In particular,
%       f(:,i) corresponds to y(:,i).
%
%   [H, f] = lrqp(A, y, Q);
%       converts the following problem (with modified norm) to QP.
%
%       minimize (1/2) * (Ax - y)' * Q * (Ax - y).
%
%       Here Q is an m x m semi-definite matrix. It can be specified
%       as either an m x m matrix, an m x 1 vector (for diagonal matrix)
%       or a scalar. Omitting Q is equivalent to setting it to be 1.
%
%   [H, f] = lrqp(A, y, Q, R);
%       converts the following ridge regression problem to QP.
%
%       minimize (1/2) (Ax - y)' * Q * (Ax - y) + (1/2) * x' * R * x
%
%       Here, R is the Tikhonov regularization coefficient matrix.
%       It is an n x n matrix, which can also be specified as an n x 1
%       vector (for diagonal matrix), or a scalar.
%

%  Created by Dahua Lin, on Nov 26, 2009
%  Modified by Dahua Lin, on Jan 2, 2010
%

%% parse and verify input arguments

if ~(isfloat(A) && ndims(A) == 2) 
    error('lr2qp:invalidarg', 'A should be a numeric matrix.');
end

if ~(isfloat(y) && ndims(y) == 2 && size(y,1) == size(A,1))
    error('lr2qp:invalidarg', ...
        'y should be a numeric matrix with same number of rows as A.');
end

%% main

if nargin < 3 || isempty(Q)
    H = A' * A;
    f = A' * y;
elseif isscalar(Q)
    H = Q * (A' * A);
    f = Q * (A' * y);
elseif size(Q,2) == 1
    H = A' * bsxfun(@times, Q, A);
    f = A' * bsxfun(@times, Q, y);
elseif ndims(Q) == 2 && size(Q,1) == size(Q,2)
    H = A' * Q * A;
    f = A' * (Q * y);
else
    error('lr2qp:invalidarg', 'The input Q is invalid.');
end

if nargin >= 4
    if isscalar(R) || isvector(R)
        H = adddiag(H, R);
    elseif ndims(R) == 2 && size(R,1) == size(R,2)
        H = H + R;
    end
end



