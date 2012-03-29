function [x, v] = hmm_decode(Pi, T, LL)
%HMM_DECODE Find the most probable state sequence for HMM
%
%   [x, v] = HMM_DECODE(Pi, T, LL);
%
%       Finds the state sequence based on an Hidden Markov Model using
%       viterbi algorithm. 
%
%       Suppose there are K states and n nodes of the chain.
%
%       Input arguments:
%       - Pi:       The initial distribution of states [K x 1]
%       - T:        The transition probability matrix [K x K]
%       - LL:       The matrix of log-likelihood values [K x n]
%
%       Output arguments:
%       - x:        The solved sequence of states
%       - v:        The joint log-likelihood of the resultant sequence.
%

% Created by Feb 2, 2012
%

%% verify input arguments

if ~(isfloat(Pi) && isreal(Pi) && ~issparse(Pi) && isvector(Pi))
    error('hmm_decode:invalidarg', ...
        'Pi should be a non-sparse real vector.');
end

K = size(T, 1);
if ~(isfloat(T) && isreal(T) && ~issparse(T) && ndims(T) == 2 && ...
        K == size(T,2))
    error('hmm_decode:invalidarg', ...
        'T should be a non-sparse real square matrix.');
end

if numel(Pi) ~= K
    error('hmm_decode:invalidarg', ...
        'The sizes of Pi and T are inconsistent.');
end

if ~(isfloat(LL) && isreal(LL) && ~issparse(LL) && ndims(LL) == 2 && ...
        size(LL,1) == K)
    error('hmm_decode:invalidarg', ...
        'LL should be a non-sparse real matrix of size K x n.');
end

%% main

if ~isa(Pi, 'double'); Pi = double(Pi); end
if ~isa(T,  'double'); T  = double(T);  end
if ~isa(LL, 'double'); LL = double(LL); end

[x, v] = chain_viterbi_cimp(log(Pi), LL, log(T), []);


