function [mu, fmsg, bmsg] = chain_bp(LA, LB)
%CHAIN_BP Discrete Belief Propagation along a chain
%
%   [mu, fmsg, bmsg] = CHAIN_BP(LA, LB);
%
%       Performs belief propagation along a chain based on the following
%       Markov model
%
%       p(x) = (1/Z) * 
%              prod_{i=1}^n a_i(x_i) * 
%              prod_{i=1}^{n-1} b(x_i, x_{i+1})
%
%       Here, we assume the second-order potentials are time-homogeneous
%
%       Suppose there are K distinct states and n nodes of the chain.
%
%       Input arguments:
%       - LA:       The logarithm of first-order potentials [K x n]
%
%       - LB:       The logarithm of second-order potential matrix,
%
%                   Note that LB can be either a K x K matrix, or a
%                   pair in form of [lb0, lb1], where lb0 is the log-
%                   potential when x_i and x_{i+1} are the same, and
%                   lb1 is the log-potential when x_i is not equal to
%                   x_{i-1}.
%
%       Output arguments:
%
%       - mu:       The marginal distribution of states [K x n]
%       - fmsg:     The (normalized) forward messages [K x (n-1)]
%       - bmsg:     The (normalized) backward messages [K x (n-1)].
%                   

% Created by Dahua Lin, on Feb 1, 2012
%

%% verify input arguments

if ~(isfloat(LA) && isreal(LA) && ~issparse(LA) && ndims(LA) == 2)
    error('chain_bp:invalidarg', 'LA should be a non-sparse real matrix.');
end
K = size(LA, 1);

if ~(isfloat(LB) && isreal(LB) && ~issparse(LB) && ...
        (isequal(size(LB), [K K]) || numel(LB) == 2))
    error('chain_bp:invalidarg', ...
        'LB should be a non-sparse real matrix of size K x K or 1 x 2.');
end

%% main

if ~isa(LA, 'double'); LA = double(LA); end
if ~isa(LB, 'double'); LB = double(LB); end

if size(LA, 2) == 1
    mu = nrmexp(LA, 1);
    fmsg = zeros(K, 0);
    bmsg = zeros(K, 0);
else
    B = exp(LB - max(LB(:)));    
    [mu, fmsg, bmsg] = chain_bp_cimp(LA, B);
end


