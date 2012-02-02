function [H, N0] = hist2d(x, y, xb, yb)
% Compute histogram of 2D data
%
%   [H, N0] = hist2d(x, y, xb, yb);
%       computes the histogram of the 2D data given by x and y, which 
%       are vectors of the same size.
%
%       The histogram bins are divided according to the xb and yb.
%       Suppose xb and yb respectively have (n+1) and (m+1) values,
%       then it divides the 2D region into n columns and m rows.
%
%       In the output, H is a matrix of size m x n, and H(i, j) is
%       the number of points in the bin at i-th row and j-th column.
%
%       N0 is the number of points falling out of the region.
%
%   H = hist2d(x, y, n, m);
%       This statement is equivalent to
%
%       hist2d(x, y, linspace(min(x), max(x), n+1), 
%                    linspace(min(y), max(y), m+1));
%

% Created by Dahua Lin, on Nov 19, 2010
%

%% verify input

if ~(isfloat(x) && isvector(x) && isreal(x) && ~issparse(x))
    error('hist2d:invalidarg', 'x should be a non-sparse real vector.');
end

if ~(isfloat(y) && isvector(y) && isreal(y) && ~issparse(y))
    error('hist2d:invalidarg', 'x should be a non-sparse real vector.');
end

if ~isequal(size(x), size(y))
    error('hist2d:invalidarg', 'x and y should have the same size.');
end

if ~(isfloat(xb) && isvector(xb))
    error('hist2d:invalidarg', 'xb should be a numeric vector.');
end

if ~(isfloat(yb) && isvector(yb))
    error('hist2d:invalidarg', 'xb should be a numeric scalar or vector.');
end

if isscalar(xb) 
    xb = linspace(min(x), max(x), xb+1);
    yb = linspace(min(y), max(y), yb+1);
end

n = numel(xb) - 1;
m = numel(yb) - 1;

%% main

ix = pieces(x, xb);
iy = pieces(y, yb);

sr = ix >= 1 & ix <= n & iy >= 1 & iy <= m;
ix = ix(sr);
iy = iy(sr);

i = sub2ind([m n], iy, ix);
H = intcount(m*n, i);

H = reshape(H, [m, n]);
N0 = numel(x) - numel(ix);


