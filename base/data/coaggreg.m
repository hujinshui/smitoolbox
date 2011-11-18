function R = coaggreg(fun, X, siz, varargin)
% Perform aggregation with multiple indices
%
%   R = coaggreg(fun, X, [K1, K2, ...], I1, I2, ...);
%       performs the aggregation specified by fun on data X with multiple
%       indices I1, I2, ...
%
%       Here, fun can be either of the following strings:
%       'sum', 'min', 'max', 'mean', 'var', 'std'.
%
%       In the input, I1, I2, ... should have the same size of X.
%
%       The returned array R will be an array of size K1 x K2 x ...
%       Concretely, R(k1, k2, ...) corresponds to the aggregation of
%       the values in X(I1==k1 & I2 == k2 & ...)
%
%       If siz is a scalar K, then R will be a column vector of size
%       K x 1.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 11, 2010
%

%% prepare

if ~(isnumeric(siz) && ndims(siz) == 2 && size(siz,1) == 1)
    error('coaggreg:invalidarg', ...
        'The 3rd argument siz should be a numeric row vector.');
end

d = numel(siz);

if numel(varargin) ~= d
    error('coaggreg:invalidarg', ...
        'The number of index arrays is invalid.');
end

I = sub2ind(siz, varargin{:});

if ~isequal(size(I), size(X))
    error('coaggreg:invalidarg', ...
        'The size of the indices arrays is invalid.');
end

%% main

if size(X, 1) ~= numel(X)
    X = X(:);
    I = I(:);
end

K = prod(siz);
R = aggreg(X, K, I, fun);

if d > 1
    R = reshape(R, siz);
end

