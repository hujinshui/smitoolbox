function dists = pwnrmdist(X1, X2, op)
%PWNRMDIST Computes pairwise normalized distances
%
%   dists = pwnrmdist(X1, X2);
%       computes pairwise normalized distances between the columns in X1
%       and X2.
%
%       Normalized distance is defined as 
%
%           d(x, y) = sqrt(2 * ( 1 - x' * y / (||x|| * ||y||) ) )
%
%       It can be proved that this distance equals the Euclidean distance
%       between normalized vectors:
%
%           d(x, y) = || (x / ||x||) - (y / ||y||)  ||  
%
%       Suppose X1 and X2 are respectively d x n1 and d x n2 matrices, then
%       dists should be an n1 x n2 matrix.
%
%       X2 can be input as an empty matrix, in which it is considered that
%       X2 = X1.
%
%   dists = pwnrmdist(X1, X2, 'square');
%       computes the squares of the pairwise normalized distances, which is 
%       just two times the the difference between one and the normalized
%       correlation.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

if nargin < 2 
    X2 = [];    
end

if isempty(X2)
    assert(isnumeric(X1) && ndims(X1) == 2, ...
        'pwnrmdist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwnrmdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pwnrmdist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end

if nargin < 3
    sqr = false;
else
    assert(strcmp(op, 'square'), ...
        'pwnrmdist:invalidarg', ...
        'The 3rd argument to nrmdist can only be ''square'' if specified.');
    
    sqr = true;
end

%% main

dists = 2 * (1 - pwnrmdot(X1, X2));

if ~sqr
    dists = sqrt(dists);
end
    


