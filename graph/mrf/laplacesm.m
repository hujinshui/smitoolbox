function X = laplacesm(Wxx, Y, wy)
% Performs Laplacian smooth based on a Gaussian MRF
%
%   The problem formulation is to minimize the following objective
%   function with respect to X:
%
%         \sum_i \sum_j Wxx(i, j) ||X(:,i) - X(:,j)||^2
%       + \sum_i wy(i) ||X(:,i) - Y(:,i)||^2.
%
%   X = laplacesm(Wxx, Y, wy);
%       solves the problems given above. Suppose X is a d x n matrix,
%       then Y should also be a d x n matrix. 
%
%       Wxx is an n x n matrix, and wy is a row vector of size 1 x n, or
%       a scalar (which indicates all wy(i) are equal).
%

% Created by Dahua Lin, on Apr 17, 2010
%

%% verify input argument

if ~(isfloat(Wxx) && isreal(Wxx) && ndims(Wxx) == 2 && size(Wxx,1) == size(Wxx,2))
    error('laplacesm:invalidarg', ...
        'Wxx should be a real symmetric matrix.');
end

n = size(Wxx, 1);

if ~(isfloat(Y) && ndims(Y) == 2 && size(Y, 2) == n)
    error('laplacesm:invalidarg', ...
        'Y should be a real matrix with n columns.');
end

if ~(isfloat(wy) && isreal(wy) && (isscalar(wy) || isequal(size(wy), [1 n])))
    error('laplacesm:invalidarg', ...
        'wy should be a real scalar or a real row vector of size 1 x n.');
end

%% main

L = laplacemat(Wxx, wy);

if isequal(wy, 1)
    YW = Y;
elseif isscalar(wy) || size(Y, 1) == 1
    YW = Y .* wy;
else
    YW = bsxfun(@times, Y, wy);
end

X = YW / L;

