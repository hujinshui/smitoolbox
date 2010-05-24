function V = pwddcrossentropy(P, Q)
%PWDDCROSSENTROPY Computes pairwise cross-entropy between distributions
%
%   V = pwddcrossentropy(P, Q);
%       computes pairwise cross-entropy values between the discrete
%       distributions in P and Q.
%
%       For two distributions with probability mass functions p and q,
%       their cross entropy is defined as
%
%           v = - sum_i p(i) * log q(i);
%
%       Let P and Q be d x m and d x n matrices respectively, then 
%       V will be an m x n matrix, with V(i, j) being the cross entropy
%       between the probabilities with their probability mass function
%       given in P(:, i) and Q(:, j).
%
%   V = pwddcrossentropy(P);
%       computes the pairwise cross-entropy between the discrete
%       distributions in columns of P. It is equivalent to 
%       pwddcrossentropy(P, P);           
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'pwddcrossentropy:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

if nargin < 2
    Q = P;
else
    assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
        'pwddcrossentropy:invalidarg', ...
        'Q should be a probability mass matrix with all elements within [0, 1]');
end

assert(size(P,1) == size(Q,1), ...
    'pwddcrossentropy:invalidarg', ...
    'P and Q should have the same number of rows.');


%% main

V = measure_mtimes_imp(P, -log(Q));
        