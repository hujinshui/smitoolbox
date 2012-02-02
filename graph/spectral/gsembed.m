function X = gsembed(g, w, d, k0)
% Spectral embedding of a graph
%
%   X = gsembed(g, w, d);
%       computes the d-dimensional embedded coordinates of a graph,
%       which correspond to the d smallest eigenvalues of the Laplacian
%       matrix.
%
%       w is the edge weights, which can be either a scalar or a 
%       vector of length m.
%
%       The output X is an n x d matrix, where X(i,:) gives the coordinate
%       of the i-th vertex.
%
%   X = gsembed(g, w, d, k0);
%       extracts the eigenvectors corresponding to (d + k0) smallest one,
%       and discards the k0 smallest, returning the remaining. 
%

% Created by Dahua Lin, on Nov 5, 2010
%

%% main

if nargin < 4
    k0 = 0;
end

L = laplacemat(g, w, 1e-10);
[X, D] = eigs(L, d+k0, 'sm'); 

dvs = diag(D);
[~, si] = sort(dvs, 1, 'ascend');
X = X(:, si(k0+1:k0+d));

