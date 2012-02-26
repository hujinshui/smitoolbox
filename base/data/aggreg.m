function R = aggreg(X, K, I, fun)
%AGGREG Index-based aggregation
%
%   R = AGGREG(X, K, I, 'sum');
%   R = AGGREG(X, K, I, 'mean');
%   R = AGGREG(X, K, I, 'min');
%   R = AGGREG(X, K, I, 'max');
%   R = AGGREG(X, K, I, 'var');
%   R = AGGREG(X, K, I, 'std');
%
%       performs aggregation of the rows or columns in X based on the 
%       indices given by I. K is the number of distinct indices.       
%
%       Suppose X is a matrix of size m x n, then I can be a vector
%       of size m x 1 or size 1 x n. 
%
%       (1) size(I) == [m, 1], then size(R) will be [K, n]:
%           R(k, :) is the aggregation of the rows in X whose
%           corresponding index in I is k. 
%
%       (2) size(I) == [1, n], then size(R) will be [m, K]
%           R(:, k) is the aggregation of the columns in X whose
%           corresponding index in I is k.
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
%       - Modified by Dahua Lin, on Mar 27, 2011
%           - re-implement the mex using bcslib
%           - now additionally supports int32 and bool.
%

%% verify input

if ~((isnumeric(X) || islogical(X)) && ~issparse(X) && ndims(X) == 2)
    error('aggreg:invalidarg', ...
        'X should be a non-sparse numeric or logical matrix.');
end
[m, n] = size(X);

if ~(isnumeric(K) && isscalar(K) && K >= 1)
    error('aggreg:invalidarg', ...
        'K should be a positive integer scalar.');
end
K = double(K);

if ~(isnumeric(I) && ~issparse(I) && isreal(I) && isvector(I))
    error('aggreg:invalidarg', ...
        'I should be a real non-sparse real vector.');
end
[mI, nI] = size(I);

if mI == m && nI == 1       % m x 1
    dim = 1;
elseif mI == 1 && nI == n   % 1 x n
    dim = 2;
else
    error('aggreg:invalidarg', ...
        'I should be a vector of size 1 x n or size m x 1.');
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
        R = ag_sum(X, K, I, dim);
        
    case 'mean'
        R = ag_mean(X, K, I, dim);        
        
    case 'min'
        R = ag_min(X, K, I, dim);
        
    case 'max'
        R = ag_max(X, K, I, dim);
        
    case 'var'
        R = ag_var(X, K, I, dim);
        
    case 'std'
        R = sqrt(ag_var(X, K, I, dim));
        
    otherwise
        error('aggreg:invalidarg', 'The fun name is invalid.');
        
end

%% sub functions

function R = ag_sum(X, K, I, dim)

R = aggreg_cimp(X, K, int32(I)-1, dim, 1); 

function R = ag_mean(X, K, I, dim)

S = ag_sum(X, K, I, dim);
if ~isfloat(S)
    S = double(S);
end
C = intcount(K, I);
R = make_mean(S, C, dim);


function R = ag_min(X, K, I, dim)

R = aggreg_cimp(X, K, int32(I)-1, dim, 2); 


function R = ag_max(X, K, I, dim)

R = aggreg_cimp(X, K, int32(I)-1, dim, 3);


function R = ag_var(X, K, I, dim)

S1 = ag_sum(X, K, I, dim);
S2 = ag_sum(X.^2, K, I, dim);
C = intcount(K, I);

if ~isfloat(S1); S1 = double(S1); end
if ~isfloat(S2); S2 = double(S2); end

E1 = make_mean(S1, C, dim);
E2 = make_mean(S2, C, dim);
R = E2 - E1.^2;
R(R < 0) = 0;



%% auxiliary function

function R = make_mean(S, C, dim)

if dim == 1
    C = C.';
    
    if size(S, 2) == 1
        R = S ./ C;
    else
        R = bsxfun(@times, S, 1./C);
    end
    
else  % dim == 2
    
    if size(S, 1) == 1
        R = S ./ C;
    else
        R = bsxfun(@times, S, 1./C);
    end
end
    

