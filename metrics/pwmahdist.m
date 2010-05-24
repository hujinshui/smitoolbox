function dists = pwmahdist(X1, X2, A, op)
%PWMAHDIST Computes the pairwise Mahalanobis distances
%
%   dists = pwmahdist(X1, X2, A);
%       computes the pairwise Mahalanobis distances between the columns in
%       X1 and X2.
%
%       Mahalanobis distances are defined as follows
%
%           d(x, y) = sqrt((x - y)' * A * (x - y));
%
%       In a d-dimensional space, both X1 and X2 should be respectively 
%       d x n1 and d x n2 matrices, and A should be a d x d symmetric matrix.
%
%       In the output, dists will be a n1 x n2 matrix, with 
%           
%           dists(i, j) = dist(X1(:,i), X2(:,j));
%
%       X2 can be input as an empty matrix, which indicates X2 = X1.
%
%   dists = pwmahdist(X1, X2, A, 'square');
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

%% parse and verify input arguments

if isempty(X2)
    assert(isnumeric(X1) && ndims(X1) == 2, ...
        'pwmahdist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwmahdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pwmahdist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end

d = size(X1, 1);
assert(ndims(A) == 2 && size(A,1) == d && size(A,2) == d, ...
    'pwmahdist:invalidarg', ...
    'A should be a d x d matrix.');

if nargin < 4
    sqr = false;
else
    assert(strcmp(op, 'square'), ...
        'pwmahdist:invalidarg', ...
        'The 4th argument to mahdist can only be ''square'' if specified.');
    
    sqr = true;
end

%% main

if isempty(X2)
    X2 = X1;
end

AX2 = A * X2;
D = (-2) * (X1' * AX2);

s1 = sum(X1 .* (A * X1), 1);
s2 = sum(X2 .* AX2, 1);

D = bsxfun(@plus, D, s1');
D = bsxfun(@plus, D, s2);

D(D < 0) = 0;

if sqr
    dists = D;
else
    dists = sqrt(D);
end


