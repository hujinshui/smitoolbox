function v = ddentropy(P)
%Computes the entropy of a discrete distribution
%
%   v = ddentropy(P);
%       computes the entropy of the discrete distributions whose probability
%       mass functions are given by P.
%
%       Given a probability mass function p, the entropy is defined as
%
%           entropy = - sum_i p(i) * log p(i)
%
%       P can be a column vector to represent the probability mass
%       function, then the function returns the entropy of the
%       corresponding distribution.       
%
%       P can be a d x n matrix to capture multiple distributions, each 
%       column represents a probability mass function over d discrete 
%       elements. In the output v is a 1 x n vector, with v(i) being the 
%       entropy of the distribution given by P(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isnumeric(P) && ndims(P) == 2 && all(P(:) >= 0 & P(:) <= 1), ...
    'ddentropy:invalidarg', ...
    'P should be a probability mass vector with all elements within [0, 1]');

%% main


if all(P(:) > 0)
    T = P .* log(P);
else
    T = zeros(size(P), class(P));
    b = P > 0;
    Pb = P(b);
    T(b) = Pb .* log(Pb);
end

v = - sum(T, 1);


