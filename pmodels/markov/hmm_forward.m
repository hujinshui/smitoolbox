function [Alpha, Lc] = hmm_forward(Pi, T, LL)
%HMM_FORWARD HMM Forward recursion
%
%   [Alpha, Lp] = HMM_FORWARD(Pi, T, LL);
%
%       Performs forward recursion to evaluate the alpha values for
%       Hidden Markov Model based inference.
%
%       Input arguments:
%       - Pi:       The initial distribution of states [K x 1]
%       - T:        The transition probability matrix [K x K]
%       - LL:       The log-likelihood matrix [K x n]
%
%       Output arguments:
%       - Alpha:    The matrix of normalized alpha-values [K x n]
%                   Alpha(k, i) := pr(z_i = k | x_1, ..., x_i)
%
%       - Lc:       The vector of log conditional probabilities [1 x n]
%                   Lc(i) = log p(x_i | x_1, ..., x_{i-1})
%
%

% Created by Dahua Lin, on Feb 1, 2012
%

%% verify input arguments

if ~(isfloat(Pi) && isreal(Pi) && ~issparse(Pi) && isvector(Pi))
    error('hmm_forward:invalidarg', ...
        'Pi should be a non-sparse real vector.');
end

K = size(T, 1);
if ~(isfloat(T) && isreal(T) && ~issparse(T) && ndims(T) == 2 && ...
        K == size(T,2))
    error('hmm_forward:invalidarg', ...
        'T should be a non-sparse real square matrix.');
end

if numel(Pi) ~= K
    error('hmm_forward:invalidarg', ...
        'The sizes of Pi and T are inconsistent.');
end

if ~(isfloat(LL) && isreal(LL) && ~issparse(LL) && ndims(LL) == 2 && ...
        size(LL,1) == K)
    error('hmm_forward:invalidarg', ...
        'LL should be a non-sparse real matrix of size K x n.');
end

%% main

if ~isa(Pi, 'double'); Pi = double(Pi); end
if ~isa(T,  'double'); T  = double(T);  end
if ~isa(LL, 'double'); LL = double(LL); end

[Alpha, Lc] = hmm_forward_cimp(Pi, T, LL);


