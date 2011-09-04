function X = gsembed(g, d, k0)
% Spectral embedding of a graph
%
%   X = gsembed(g, d);
%       computes the d-dimensional embedded coordinates of a graph,
%       which correspond to the d smallest eigenvalues of the Laplacian
%       matrix.
%
%       The output X is an n x d matrix, where X(i,:) gives the coordinate
%       of the i-th vertex.
%
%   X = gsembed(g, d, k0);
%       extracts the eigenvectors corresponding to (d + k0) smallest one,
%       and discards the k0 smallest, returning the remaining. 
%

% Created by Dahua Lin, on Nov 5, 2010
%

%% main

if nargin < 3
    k0 = 0;
end

L = laplacemat(g, 1e-10);
[X, D] = eigs(L, d+k0, 'sm'); %#ok<NASGU>

if k0 > 0
    X = X(:, 1:d);
end

