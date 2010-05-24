function dists = mahdist(X1, X2, A, op)
%MAHDIST Computes the Mahalanobis distances
%
%   dists = mahdist(X1, X2, A);
%       computes the Mahalanobis distances between the corresponding
%       columns in X1 and X2.
%
%       Mahalanobis distances are defined as follows
%
%           d(x, y) = sqrt((x - y)' * A * (x - y));
%
%       In a d-dimensional space, both X1 and X2 should be d x n matrices,
%       and A should be a d x d symmetric matrix.
%
%       In the output, dists will be a 1 x n vector, with 
%           
%           dists(i) = dist(X1(:,i), X2(:,i));
%
%   dists = mahdist(X1, X2, A, 'square');
%       computes the squares of the Mahalanobis distances, i.e. not
%       taking the square root in the above formula.
%
%   Remarks
%       - It is the caller's responsibility to ensure A is symmetric, the
%         function in itself does not perform the checking for efficiency.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% main

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'mahdist:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'mahdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

assert(ndims(A) == 2 && size(A,1) == d && size(A,2) == d, ...
    'mahdist:invalidarg', ...
    'A should be a d x d matrix.');

if nargin < 4
    sqr = false;
else
    assert(strcmp(op, 'square'), ...
        'mahdist:invalidarg', ...
        'The 4th argument to mahdist can only be ''square'' if specified.');
    
    sqr = true;
end


%% main

D = X1 - X2;
dists = sum(D .* (A * D), 1);
dists(dists < 0) = 0;

if ~sqr
    dists = sqrt(dists);
end
