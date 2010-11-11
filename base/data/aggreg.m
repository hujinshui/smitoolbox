function R = aggreg(X, K, I, fun)
% Perform index-based aggregation
%
%   R = aggreg(X, K, I, 'sum');
%   R = aggreg(X, K, I, 'mean');
%   R = aggreg(X, K, I, 'min');
%   R = aggreg(X, K, I, 'max');
%   R = aggreg(X, K, I, 'var');
%   R = aggreg(X, K, I, 'std');
%
%       performs aggregation of the rows or columns in X based on the 
%       indices given by I. K is the number of distinct indices.       
%
%       Suppose X is a matrix of size m x n, then I can be a vector
%       of size m x 1 or size 1 x n. 
%       If size(I) == [m, 1], then the output matrix R is of size K x n,
%       such that R(k, :) is the aggregation of the rows in X whose
%       corresponding index in I is k. 
%       If size(I) == [1, n], then the output matrix R is of size m x K,
%       such that R(:, k) is the aggregation of the columns in X whose
%       corresponding index in I is k.
%
%       If the 4th argument is omitted, it is set to 'sum' by default.
%
%   Remarks
%   -------
%       - For the index k which corresponds to no elements in X, the
%         corresponding output value will be set to 0, NaN, Inf, and
%         -Inf, respectively for sum, mean, max, and min.
%
%       - In computing the variance, we use n rather than n-1 to normalize
%         the output. In other words, it results in the maximum likelihood
%         estimation, instead of the unbiased estimation.
%
%       - This function is partly similar to the built-in function
%         accumarray. There are two main differences: 
%           (1) it is based on optimized C++ code, and thus runs much
%               faster for the aggregation functions that it can support
%           (2) it supports aggregation of rows and columns, where
%               accumarray only supports aggregation of scalars.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 11, 2010
%

%% verify input

if ~(isnumeric(X) && ~issparse(X) && isreal(X) && ndims(X) == 2)
    error('aggreg:invalidarg', ...
        'X should be a non-sparse real matrix.');
end
[m, n] = size(X);

if ~(isnumeric(K) && isscalar(K) && K >= 1)
    error('aggreg:invalidarg', ...
        'K should be a positive integer scalar.');
end
K = double(K);

if ~(isnumeric(I) && ~issparse(I) && isreal(I) && ...
    (isequal(size(I), [m 1]) || isequal(size(I), [1 n])) )
    error('aggreg:invalidarg', ...
        'I should be a non-sparse real vector of size m x 1 or 1 x n.');
end

if nargin < 4
    fun = 'sum';
else
    if ~ischar(fun)
        error('aggreg:invalidarg', 'The fun name must be a string.');
    end
end

%% main

switch fun
    case 'sum'
        R = ag_sum(X, K, I);
        
    case 'mean'
        R = ag_mean(X, K, I);        
        
    case 'min'
        R = ag_min(X, K, I);
        
    case 'max'
        R = ag_max(X, K, I);
        
    case 'var'
        R = ag_var(X, K, I);
        
    case 'std'
        R = sqrt(ag_var(X, K, I));
        
end


%% sub functions

function R = ag_sum(X, K, I)

R = aggreg_cimp(X, K, int32(I)-1, 1); 

function R = ag_mean(X, K, I)

S = ag_sum(X, K, I);
C = intcount(K, I);
R = make_mean(S, C, I);


function R = ag_min(X, K, I)

R = aggreg_cimp(X, K, int32(I)-1, 2); 


function R = ag_max(X, K, I)

R = aggreg_cimp(X, K, int32(I)-1, 3);


function R = ag_var(X, K, I)

S1 = ag_sum(X, K, I);
S2 = ag_sum(X.^2, K, I);
C = intcount(K, I);

E1 = make_mean(S1, C, I);
E2 = make_mean(S2, C, I);
R = E2 - E1.^2;
R(R < 0) = 0;



%% auxiliary function

function R = make_mean(S, C, I)

if size(I, 1) ~= 1
    C = C.';
    if size(S, 2) == 1
        R = S ./ C;
    else
        R = bsxfun(@times, S, 1./C);
    end
else
    if size(S, 1) == 1
        R = S ./ C;
    else
        R = bsxfun(@times, S, 1./C);
    end
end


