function D = ddjsdiv(P, Q)
%DDJSDIV Computes the J-S Divergence of discrete distributions
%
%   D = ddjsdiv(P, Q);
%       computes the Jensen-Shannon divergence between discrete
%       distributions, whose probability mass functions are given in 
%       P and Q.
%
%       Jensen-Shannon divergence, also known as Jeffrey divergence,
%       is a "symmetric version" of the well-known Kullback-Leibler 
%       divergence. J-S divergence between p and q is defined as
%
%       JSD(p || q) = (1/2) sum_i p(i) log p(i)/u(i) + q(i) log q(i)/u(i)
%
%       Here, u = (p + q)/2, which is the "average distribution".
%
%       Compared to K-L divergence, J-S divergence is symmetric and
%       numerically stable.
%
%       P and Q should be matrices of the same size, let it be m x n.
%       Then D will be a 1 x n row vector, with D(i) giving the J-S
%       divergence between the distributions characterized by the
%       probability mass function P(:, i) and Q(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'ddjsdiv:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

assert(isfloat(Q) && ndims(Q) == 2 && all(Q(:) >= 0 & Q(:) <= 1), ...
    'ddjsdiv:invalidarg', ...
    'Q should be a probability mass matrix with all elements within [0, 1]');

assert(size(P,1) == size(Q,1) && size(P,2) == size(Q,2), ...
    'ddjsdiv:invalidarg', ...
    'P and Q should have the same size.');
       
%% main

D = ddentropy(0.5 * (P + Q)) - 0.5 * (ddentropy(P) + ddentropy(Q));

