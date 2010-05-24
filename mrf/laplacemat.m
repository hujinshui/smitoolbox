function L = laplacemat(W, a)
% Compute the (regularized) Laplacian matrix of a graph
%
%   L = laplacemat(W);
%   L = laplacemat(W, a);
%       computes the Laplacian matrix defined as follows:
%
%           L = diag(sum(W) + a) - Wxx.
%
%       Here, W the affinity matrix of size n x n. And a can be 
%       in either of the following form:
%       - a vector of length n
%       - a scalar
%       - omited or empty, when a is 0.
%
%       The Laplacian matrix plays an important role in graph-based
%       optimization.
%

% Created by Dahua Lin, on Apr 17, 2010
%

%% verify input arguments

if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,1) == size(W,2))
    error('laplacemat:invalidarg', 'W should be a real symmetric matrix.');
end

n = size(W, 1);

if nargin < 2 || isempty(a)
    a = 0;
else
    if ~(isfloat(a) && (isscalar(a) || (isvector(a) && numel(a) == n)))
        error('laplacemat:invalidarg', ...
            'a should be a real scalar or a real vector of length n.');
    end
    
    if size(a, 2) > 1   % turn a into a column vector
        a = a.';
    end
end

%% main

dv = sum(W, 2);

if issparse(dv)
    dv = full(dv);
end

if ~isequal(a, 0)
    dv = dv + a;
end

if issparse(W)
    L = spdiag(dv) - W;
else
    di = 1 : (0:n-1) * (n+1);
    L = -W;
    L(di) = L(di) + a;
end

