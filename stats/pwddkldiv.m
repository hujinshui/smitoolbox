function D = pwddkldiv(P, Q)
%PWDDKLDIV Computes pairwise K-L divergences
%
%   D = pwddkldiv(P, Q);
%       computes pairwise Kullback-Leibler divergence between the discrete
%       distributions in P and Q.
%
%       K-L divergence between two probability distributions p and q is 
%       defined as
%
%           KL(p || q) = sum_i p(i) * log ( p(i) / q(i) )
%
%       The K-L divergence has the following identity:
%
%           KL(p || q) = H(p, q) - H(p)
%
%       here, H(p, q) is the cross-entropy between p and q, H(p) is p's 
%       entropy.
%
%       Let P and Q be d x m and d x n matrices respectively, then 
%       D will be an m x n matrix, with D(i, j) being the K-L divergence
%       between the probabilities with their probability mass function
%       given in P(:, i) and Q(:, j).
%
%   D = pwddkldiv(P);
%       computes the pairwise K-L divergence between the discrete
%       distributions in columns of P. It is equivalent to 
%       pwddkldiv(P, P);           
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'pwddkldiv:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

if nargin < 2
    Q = P;
else
    assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
        'pwddkldiv:invalidarg', ...
        'Q should be a probability mass matrix with all elements within [0, 1]');
end

assert(size(P,1) == size(Q,1), ...
    'pwddkldiv:invalidarg', ...
    'P and Q should have the same number of rows.');

%% main

D = bsxfun(@minus, pwddcrossentropy(P, Q), ddentropy(P)');

