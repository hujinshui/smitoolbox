function r = sum_xlogy(X, Y, dim)
% Compute the sum of the terms x * log(y) in a safe way
%
%   r = sum_xlogy(X, Y);
%   
%       Here, X and Y should be matrices or vectors of the same type.
%
%       If X and Y are both vectors, then r is a scalar.
%       Otherwise, if they are m x n matrices, then the sum is taken 
%       along the first dimension, resuling in r of size 1 x n.
%
%       Note that in this function, 0 * log(y) is zero for every y,
%       in particular, 0 * log(0) is regarded as 0. This is very useful
%       for computing information theoretical quantities.
%
%   r = sum_xlogy(X, Y, dim);
%
%       Perform the sum of x * log(y) along the specified dimension.
%       Here, dim should be either 1 or 2.
%
%       Note that the sizes of X and Y along the specified dimension
%       should be the same. However, the following input is allowed:
%       one is vector, while the other is a matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on Mar 28, 2011
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2 && isfloat(Y) && isequal(size(X), size(Y)))
    error('sum_xlogy:invalidarg', ...
        'X and Y should be numeric matrices of the same size.');
end

%% main

Y(X == 0) = 1;

if nargin < 3
    r = sum(X .* log(Y));
else
    r = sum(X .* log(Y), dim);
end

