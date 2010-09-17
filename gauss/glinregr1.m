function [c1a, c2a] = glinregr1(X, y, w)
% Get posterior update for linear regression with scalar response
%
%   Consider the likelihood model formalized as follows:
%
%       y ~ N(a' * x, sigma^2)
%    
%   This function computes the posterior update to the 1st order and
%   2nd order coefficients of the Gaussian prior of a.
%
%   [c1a, c2a] = glinregr1(X, y);
%   [c1a, c2a] = glinregr1(X, y, w);
%       captures the observations in X and y. Here, X(:,i) is the i-th
%       input vector, and y(i) is the corresponding response. 
%       w is the sample weight, which equals 1 / sigma^2.
%       
%       It outputs c1a and c2a, which are used to update coef1 and coef2
%       respectively for a Gaussian prior distribution.
%
%       Formally, the formulas to compute c1a and c2a are given as
%       follows:
%
%           c1a = sum_i w_i * (y_i * x_i);
%           c2a = sum_i w_i * (x_i * x_i');
%
%       Suppose X is a matrix of size d x n, then y should be a row vector
%       of size 1 x n. w can be a scalar (if all samples share the same
%       weight) or a row vector of size 1 x n. If w is omitted, all samples
%       have the same weight 1.

% Created by Dahua Lin, on Sep 16, 2010
%


%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('glinrepr1:invalidarg', 'X should be a numeric matrix.');
end

n = size(X, 2);

if ~(isfloat(y) && ndims(y) == 2 && size(y,1) == 1 && size(y,2) == n)
    error('glinrepr1:invalidarg', 'y should be a numeric row vector of size 1 x n.');
end

if nargin < 3
    w = 1;
else
    if ~(isfloat(w) && (isscalar(w) || (ndims(w) == 2 && size(w,1) == 1 && size(w,2) == n)))
        error('glinrepr1:invalidarg', ...
            'w should be either a scalar or a numeric row vector of size 1 x n.');
    end
end


%% main

if isscalar(w)
    
    c1a = X * y';
    c2a = X * X';
    
    if w ~= 1
        c1a = w * c1a;
        c2a = w * c2a;
    end
    
else
    
    c1a = X * (w .* y)';
    c2a = X * bsxfun(@times, X, w)';
    
    c2a = 0.5 * (c2a + c2a');
    
end

c2a = gsymat(c2a);

