function R = aggreg_percol(X, K, I, fun)
%AGGREG_PERCOL Per-column index-based aggregation
%
%   R = AGGREG_PERCOL(X, K, I, 'sum');
%   R = AGGREG_PERCOL(X, K, I, 'mean');
%   R = AGGREG_PERCOL(X, K, I, 'min');
%   R = AGGREG_PERCOL(X, K, I, 'max');
%   R = AGGREG_PERCOL(X, K, I, 'var');
%   R = AGGREG_PERCOL(X, K, I, 'std');
%
%       performs per-column aggregation of X based on the indices given
%       by I. K is the number of distinct indices.
%
%       X and I should be matrices of the same size. Suppose this size
%       is m x n, then R will be a matrix of size K x n, such that
%
%           R(k, i) = afun(X(I(:,i) == k, i));
%
%       Here, afun is the given aggregation function.
%       
%       If the 4th argument is omitted, it is set to 'sum' by default.
%

% Created by Dahua Lin, on Feb 26, 2012
%

%% verify inputs

if ~((isnumeric(X) || islogical(X)) && ~issparse(X) && ismatrix(X))
    error('aggreg_percol:invalidarg', ...
        'X should be a non-sparse numeric or logical matrix.');
end

if ~(isnumeric(K) && isscalar(K) && K >= 1)
    error('aggreg_percol:invalidarg', ...
        'K should be a positive integer scalar.');
end
K = double(K);

if ~(isnumeric(I) && ~issparse(I) && isreal(I) && ismatrix(I))
    error('aggreg_percol:invalidarg', ...
        'I should be a real non-sparse real vector.');
end

if ~(size(X,1) == size(I,1) && size(X,2) == size(I,2))
    error('aggreg_percol:invalidarg', ...
        'X and I should be matrices of the same size.');
end

if nargin < 4
    fun = 'sum';
else
    if ~ischar(fun)
        error('aggreg_percol:invalidarg', 'The fun name must be a string.');
    end
end

%% main

switch fun
    case 'sum'
        R = agx_sum(X, K, I);
        
    case 'mean'
        R = agx_mean(X, K, I);        
        
    case 'min'
        R = agx_min(X, K, I);
        
    case 'max'
        R = agx_max(X, K, I);
        
    case 'var'
        R = agx_var(X, K, I);
        
    case 'std'
        R = sqrt(agx_var(X, K, I));
        
    otherwise
        error('aggreg:invalidarg', 'The fun name is invalid.');
        
end


%% Core functions

function R = agx_sum(X, K, I)

R = aggreg_percol_cimp(X, K, int32(I)-1, 1); 


function R = agx_mean(X, K, I)

S = agx_sum(X, K, I); 
if ~isfloat(S)
    S = double(S);
end
N = intcount(K, I, 1);
R = S ./ N;


function R = agx_min(X, K, I)

R = aggreg_percol_cimp(X, K, int32(I)-1, 2); 


function R = agx_max(X, K, I)

R = aggreg_percol_cimp(X, K, int32(I)-1, 3);


function R = agx_var(X, K, I)

S1 = agx_sum(X, K, I); 
S2 = agx_sum(X.^2, K, I);

if ~isfloat(S1)
    S1 = double(S1);
end
if ~isfloat(S2)
    S2 = double(S2);
end
       
N = intcount(K, I, 1);
EX = S1 ./ N;
EX2 = S2 ./ N;

R = EX2 - EX.^2;
R(R < 0) = 0;





