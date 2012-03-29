function Beta = hmm_backward(T, LL, Lc)
%HMM_BACKWARD HMM Backward recursion
%
%   Beta = HMM_BACKWARD(T, LL, Lc);
%
%       Performs forward recursion to evaluate the alpha values for
%       Hidden Markov Model based inference.
%
%       Input arguments:
%       - T:        The transition probability matrix [K x K]
%       - LL:       The log-likelihood matrix [K x n]
%       - Lc:       The log conditional likelihood vector [1 x n]
%                   (This is the 2nd output argument of hmm_forward).
%
%       Output arguments:
%       - Beta:     The vector of normalized beta values [K x n]
%
%

% Created by Dahua Lin, on Feb 1, 2012
%

%% verify input arguments

K = size(T, 1);
if ~(isfloat(T) && isreal(T) && ~issparse(T) && ndims(T) == 2 && ...
        K == size(T,2))
    error('hmm_backward:invalidarg', ...
        'T should be a non-sparse real square matrix.');
end

if ~(isfloat(LL) && isreal(LL) && ~issparse(LL) && ndims(LL) == 2 && ...
        size(LL,1) == K)
    error('hmm_backward:invalidarg', ...
        'LL should be a non-sparse real matrix of size K x n.');
end
n = size(LL, 2);

if ~(isfloat(Lc) && isreal(Lc) && isequal(size(Lc), [1 n]))
    error('hmm_backward:invalidarg', ...
        'Lc should be a non-sparse real vector of size 1 x n.');
end


%% main

if ~isa(T,  'double'); T  = double(T);  end
if ~isa(LL, 'double'); LL = double(LL); end
if ~isa(Lc, 'double'); Lc = double(Lc); end

Beta = hmm_backward_cimp(T, LL, Lc);


