function Pr = regulateprob(P, a, dim)
%REGULATEPROB Regulates the probabilities by adding prior counts
%
%   Pr = regulateprob(P, a)
%       regulates the probabilities in P by adding prior counts.
%
%       The probability mass functions with zero entries are instable and
%       difficult to deal with in many cases, like those involving the
%       computation of entropy or KL-divergence.
%
%       This function implements a simple regulation technique that
%       slightly modifies the original probability mass function to
%       make it easier to handle.
%
%       Given a probability mass function p, and a prior count a, then
%       pr, the modified function, is given by
%
%           pr(i) = max(p(i) + a, 0) / sum_i max(p(i) + a, 0)
%
%       When a is positive, it ensures that all entries of p become
%       positive.
%
%       When a is negative, it attempts to bring some small entries to
%       zeros.
%
%       In the input, P can be a vector to characterize a distribution,
%       then Pr will be a vector of the same size to characterize the
%       regulated distribution.
%
%       P can also be a matrix of size m x n with m > 1, then each column 
%       of P is considered to represent a distribution. 
%
%   Pr = regulateprob(P, a, 1)
%       regulates the probabilities in P, with each column representing a
%       distribution.
%
%   Pr = regulateprob(P, a, 2)
%       regulates the probabilities in P, with each row representing a
%       distribution.
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'regulateprob:invalidarg', ...
    'P should be a probability mass matrix with all elements within [0, 1]');

assert(isfloat(a) && isscalar(a), 'regulateprob:invalidarg', ...
    'a should be a floating-point scalar.');

if nargin < 3
    if size(P, 1) > 1
        dim = 1;
    else
        dim = 2;
    end
    
else
    assert(isscalar(dim) && (dim == 1 || dim == 2), ...
        'regulateprob:invalidarg', 'dim should be either 1 or 2.');
    
end

%% main

if a > 0
    Pr = P + a;

elseif a < 0
    Pr = max(P + a, 0);

end

Pr = bsxfun(@times, Pr, 1 ./ sum(Pr, dim));


