function [R, A] = knnc(D, M, L, K, wfun, op)
%KNNC K-nearest-neighbor classification
%
%   R = KNNC(D, M, L, K);
%   R = KNNC(D, M, L, K, wfun);
%   R = KNNC(D, M, L, K, [], 'dis');
%   R = KNNC(D, M, L, K, wfun, 'dis');
%   R = KNNC(D, M, L, K, [], 'sim');
%   R = KNNC(D, M, L, K, wfun, 'sim');
%
%       Performs K-nearest-neighbor classification. 
%
%       Input arguments:
%       - D:        The matrix of distances or similarities.
%                   
%                   Suppose there are n query samples and m reference
%                   samples, then D should have size m x n. Particularly,
%                   D(i, j) is the distance between the j-th query to
%                   the i-th reference.
%
%                   If the 6th argument is omited or set to 'dis', then
%                   D contains distances (i.e. smaller values indicate
%                   closer). If the 6th argument is set to 'sim', then
%                   D contains similarities (i.e. higher values indiate
%                   closer).
%
%       - M:        The number of distinct classes. 
%
%       - L:        The labels of reference samples, which should be
%                   a numeric vector of length m. The values in L should
%                   be integers in {1, ..., M}.
%
%       - K:        The number of neighbors for each sample.
%
%       - wfun:     The weighting function, which, if given, will be
%                   invoked as
%
%                       W = wfun(Dr);
%
%                   Here, Dr is the matrix of sorted distance or 
%                   similarity values (from the closest to the K-th 
%                   closest). The size of Dr is K x n. In the output,
%                   W should also be a matrix of size K x n.
%
%                   If wfun is omitted or input as [], then all the K
%                   neighbors have equal contribution.
%
%       Output arguments:
%       - R:        The resultant vector of class labels, of size 1 x n.
%
%
%   [R, A] = KNNC( ... );
%       
%       Additionally returns A, the class-wise accumulated scores for
%       each pair of class/sample. 
%
%       Here, A is a matrix of size K x n. A(k, i) is the number of
%       neighbors of sample i that belong to class k (if wfun is not used),
%       or the total neighbor weights associated with class k.
%

% Created by Dahua Lin, on Jan 21, 2012
%

%% verify input

if ~(isfloat(D) && isreal(D) && ndims(D) == 2)
    error('knnc:invalidarg', 'D should be a real matrix.');
end
[m, n] = size(D);
if m < 2
    error('knnc:invalidarg', 'D should contain at least two rows.');
end

if ~(isnumeric(M) && isscalar(M) && M == fix(M) && M >= 1)
    error('knnc:invalidarg', 'M should be a positive integer.');
end

if ~(isnumeric(L) && isvector(L) && numel(L) == m)
    error('knnc:invalidarg', 'L should be a vector of length m.');
end

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1 && K <= n)
    error('knnc:invalidarg', 'K should be a positive integer in [1, n].');
end

if nargin < 5 || isempty(wfun)
    use_wfun = 0;
else
    if ~isa(wfun, 'function_handle')
        error('knnc:invalidarg', 'wfun should be a function handle.');
    end
    use_wfun = 1;
end

if nargin < 6
    dir = 'min';
else
    if strcmpi(op, 'dis')
        dir = 'min';
    elseif strcmpi(op, 'sim')
        dir = 'max';
    else
        error('knnc:invalidarg', 'The 6th argument is invalid.');
    end
end
        

%% main

[Dr, I] = top_k(D, dir, K);
Lr = L(I);

if ~use_wfun
    A = intcount(M, Lr, 1);
else
    W = wfun(Dr);
    assert(isequal(size(W), [K n]));
    if ~isa(W, 'double')
        W = double(W);
    end    
    
    A = aggreg_percol(W, M, Lr, 'sum');
end

[~, R] = max(A, [], 1);



