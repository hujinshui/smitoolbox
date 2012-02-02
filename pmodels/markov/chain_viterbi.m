function [x, v] = chain_viterbi(A, B)
%CHAIN_VITERBI Viterbi algorithm on a Hidden Markov chain
%
%   [x, v] = CHAIN_VITERBI(A, B);
%
%       Performs the Viterbi algorithm along a chain that finds a sequence
%       of states that maximizes the following objective
%
%           sum_{i=1}^n a_i(x_i) + sum_{i=1}^{n-1} b(x_i, x_{i+1})
%
%
%       Suppose there are K distinct states and n nodes of the chain.
%
%       Input arguments:
%
%       - A:        The matrix of first-order (additive) potentials 
%                   [size: K x n].
%
%       - B:        The matrix of second-order (additive) potentials,
%                   which can be in either of the following two forms:
%                   - K x K matrix, then b(u, v) = B(u, v).
%                   - 1 x 2 vector, then b(u, v) = B(1) when u == v
%                     and b(u, v) = B(2) when u != v.
%
%       Output arguments:
%       - x:        The resultant sequence of states [1 x n].
%
%       - v:        The objective (total potential) of the optimal state
%                   sequence.
%

% Created by Dahua Lin, on Feb 2, 2012
%

%% verify input arguments

if ~(isfloat(A) && isreal(A) && ~issparse(A) && ndims(A) == 2)
    error('chain_viterbi:invalidarg', ...
        'A should be a non-sparse real matrix.');
end

if ~( isfloat(B) && isreal(B) && ~issparse(B) && ...
        (isequal(size(B), [K K]) || numel(B) == 2) )
    error('chain_viterbi:invalidarg', ...
        'B should be a non-sparse real matrix of size K x K or 1 x 2.');
end

%% main

if ~isa(A, 'double'); A = double(A); end
if ~isa(B, 'double'); B = double(B); end

[x, v] = chain_viterbi(A, B);



        