function P = binoprobs(n, p)
% Compute probability of binomial distribution
%
%   P = binoprobs(n, p);
%       Compute the probability mass function of binomial distribution.
%       In the input, n is a positive integer representing the total
%       number of trials, p is the probability of success of each trial.
%       
%       Then in the output, P is a row vector of size 1 x (n+1), where
%       P(i) is the probability of having i - 1 successes in n trials.
%
%       p can also be a vector with m elements that give different
%       probabilities of success. In this case, P is a matrix of size
%       m x (n+1), where P(k, :) corresponds to the probabilities derived
%       with success rate p(k).
%

% Created by Dahua Lin, on Apr 16, 2010
%

%% verify input arguments

if ~(isfloat(n) && isscalar(n) && n == fix(n) && n >= 1)
    error('binoprobs:invalidarg', ...
        'n should be a positive integer scalar.');
end

if ~(isfloat(p) && isreal(p) && isvector(p))
    error('binoprobs:invalidarg', ...
        'p should be a real vector.');
end

if size(p, 2) > 1
    p = p.';
end
m = numel(p);


%% main

fv = [1 cumprod(1:n)]; % the list of values of factorials
k = 0 : n;
fn = fv(n+1);

if m == 1
    a = (p .^ k) ./ fv;
    b = ((1 - p) .^ (n - k)) ./ fv(end:-1:1);
    P = fn * (a .* b);
else
    A = bsxfun(@times, bsxfun(@power, p, k), 1 ./ fv);
    B = bsxfun(@times, bsxfun(@power, 1-p, n-k), 1 ./ fv(end:-1:1));
    P = fn * (A .* B);
end

