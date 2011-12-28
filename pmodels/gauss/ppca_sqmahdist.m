function dists = ppca_sqmahdist(M, X)
% Compute squared Mahalanobis distances to the center of PPCA model
%
%   dists = ppca_sqmahdist(M, X);
%       computes the squared Mahalanobis distances from the samples 
%       in X to the center of the PPCA model M, which is defined by
%
%           (x - mu)' * inv(C) * (x - mu).
%
%       Suppose X has n columns (each column is sample), then dists
%       is a vector of size 1 x n.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~is_ppca(M)
    error('ppca_sqmahdist:invalidarg', 'M should be a PPCA struct.');
end

if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == M.d)
    error('ppca_sqmahdist:invalidarg', ...
        'X should be a real matrix with d rows.');
end

%% main

% center

mu = M.mu;
if ~isequal(mu, 0)
    X = bsxfun(@minus, X, mu);
end

% decompose

B = M.B;
C = B' * X;
E = X - B * C;

% evaluate

v2 = M.se^2;
v1 = M.s.^2 + v2;

dists = (1 ./ v1) * (C.^2) + (1/v2) * sum(E.^2, 1);


