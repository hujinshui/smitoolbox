function v = ddcrossentropy(P, Q)
%DDCROSSENTROPY Computes the cross-entropy of discrete distributions
%
%   v = ddcrossentropy(P, Q);
%       computes the cross entropy between discrete distributions given
%       by P and Q.
%
%       For two distributions with probability mass functions p and q,
%       their cross entropy is defined as
%
%           v = - sum_i p(i) * log q(i);
%
%       In the input, P and Q should be numeric matrices of the same 
%       size d x n.
%
%       If n == 1, i.e. P and Q are both column vectors, then they
%       respectively characterize the distributions p and q, then v 
%       will be their cross entropy.
%
%       If n > 1, then each column of P and Q characterizes a distribution,
%       then v will be a 1 x n row vector, with v(i) gives the cross
%       entropy between P(:, i) and Q(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'ddcrossentropy:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
    'ddcrossentropy:invalidarg', ...
    'Q should be a probability mass matrix with all elements within [0, 1]');

assert(size(P,1) == size(Q,1) && size(P,2) == size(Q,2), ...
    'ddcrossentropy:invalidarg', ...
    'P and Q should have the same size.');

%% main

if all(P > 0)
    T = P .* log(Q);
    
else
    T = zeros(size(P), class(P));
    b = P > 0;    
    T(b) = P(b) .* log(Q(b));
    
end

v = - sum(T, 1);
