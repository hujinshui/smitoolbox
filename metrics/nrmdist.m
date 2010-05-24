function dists = nrmdist(X1, X2, op)
%NRMDIST Computes normalized distances between corresponding vectors
%
%   dists = nrmdist(X1, X2);
%       computes normalized distances between corresponding vectors in X1
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
%       X1 and X2 should be both d x n matrices, then dists will be a 1 x n
%       vector, with dists(i) being the normalized distance between the
%       i-th vector in X1 and that in X2.
%
%   dists = nrmdist(X1, X2, 'square');
%       computes the square of normalized distances, which is just two
%       times the the difference between one and the normalized
%       correlation.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'nrmdist:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'nrmdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

if nargin < 3
    sqr = false;
else
    assert(strcmp(op, 'square'), ...
        'nrmdist:invalidarg', ...
        'The 3rd argument to nrmdist can only be ''square'' if specified.');
    
    sqr = true;
end


%% main

s1 = sum(X1 .* X1, 1);
s2 = sum(X2 .* X2, 1);

dists = 2 * (1 - sum(X1 .* X2, 1) ./ sqrt(s1 .* s2));      
dists(dists < 0) = 0;    

if ~sqr
    dists = sqrt(dists);
end


