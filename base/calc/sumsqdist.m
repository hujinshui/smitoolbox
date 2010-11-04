function s = sumsqdist(X, Y, W)
% Compute weighted sum of square distances
%
%   s = sumsqdist(X, Y, W);
%       computes the weighted sum of squared distance as follows:
%
%       s = \sum_{i=1}^m \sum_{j=1}^n W(i,j) ||X(:,i) - Y(:,j)||^2.
%
%       Here, X is a d x m matrix, Y is a d x n matrix, and 
%       W is an m x n matrix. 
%
%   s = sumsqdist(X, [], W);
%       computes the weighted sum with X = Y and a symmetric matrix W.
%       More efficient implementation will be used for this case.
%

% Created by Dahua Lin, on Apr 17, 2010
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('sumsqdist:invalidarg', 'X should be a real matrix.');
end

if isempty(Y)
    use_sym = true;
    Y = X;
else
    if ~(isfloat(Y) && isreal(Y) && ndims(Y) == 2)
        error('sumsqdist:invalidarg', 'Y should be a real matrix.');
    end
    use_sym = false;
end

if size(X, 1) ~= size(Y, 1)
    error('sumsqdist:invalidarg', ...
        'x-vectors and y-vectors have different dimension.');
end

d = size(X, 1);
m = size(X, 2);
n = size(Y, 2);

if ~(isfloat(W) && isreal(W) && isequal(size(W), [m n]))
    error('sumsqdist:invalidarg', 'W should be an m x n real matrix.');
end

%% main

if d == 1
    sx = sum((X.^2) .* sum(W, 2).');
    
    if use_sym
        sy = sx;
    else
        sy = sum((Y.^2) .* sum(W, 1));
    end
    
    sxy = sum((X * W) .* Y);
    
else
    sx = sum(sum(bsxfun(@times, X, sum(W, 2).') .* X, 1), 2);

    if use_sym
        sy = sx;
    else
        sy = sum(sum(bsxfun(@times, Y, sum(W, 1)) .* Y, 1), 2);
    end
    
    sxy = sum(sum((X * W) .* Y, 1), 2);
end

s = sx + sy - 2 * sxy;

