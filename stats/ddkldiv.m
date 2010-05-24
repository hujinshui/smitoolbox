function D = ddkldiv(P, Q)
%DDKLDIV Computes K-L divergence between discrete distributions
%
%   D = ddkldiv(P, Q)
%       computes Kullback-Leibler (K-L) divergence between distribution(s),
%       whose probability mass functions are given in P and Q.
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
%       P and Q should be matrices of the same size, let it be m x n.
%       Then D will be a 1 x n row vector, with D(i) giving the K-L
%       divergence between the distributions characterized by the
%       probability mass function P(:, i) and Q(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'ddkldiv:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
    'ddkldiv:invalidarg', ...
    'Q should be a probability mass matrix with all elements within [0, 1]');

assert(size(P,1) == size(Q,1) && size(P,2) == size(Q,2), ...
    'ddkldiv:invalidarg', ...
    'P and Q should have the same size.');

%% main

if all(P(:) > 0)
    
    T = P .* (log(P) - log(Q));
    
else
    
    T = zeros(size(P), class(P));
    b = P > 0;
    Pb = P(b);
    T(b) = Pb .* (log(Pb) - log(Q(b)));
    
end
    
D = sum(T, 1);

